#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <curl/curl.h>
#include <opencv2/opencv.hpp>

struct GpxPoint {
    double lat;
    double lon;
    double ele;
    std::string time;
};

// Forward declaration for the fallback function
cv::Mat downloadOpenStreetMapTile(double north, double south, double east, double west, int width, int height);

// Convert WGS84 (lat/lon) to Web Mercator coordinates
void latLonToWebMercator(double lat, double lon, double& x, double& y) {
    // Convert to radians
    lat = lat * M_PI / 180.0;
    lon = lon * M_PI / 180.0;
    
    // Web Mercator projection (EPSG:3857)
    x = (lon + M_PI) / (2.0 * M_PI);
    y = (1.0 - log(tan(lat) + 1.0 / cos(lat)) / M_PI) / 2.0;
}

// Convert Web Mercator coordinates to pixel coordinates
cv::Point webMercatorToPixel(double x, double y, double west, double east, double north, double south, int width, int height) {
    // Convert Web Mercator coordinates to pixel coordinates
    int pixelX = static_cast<int>(std::round(((x - west) / (east - west)) * width));
    int pixelY = static_cast<int>(std::round(((y - north) / (south - north)) * height));
    
    return cv::Point(pixelX, pixelY);
}

// Function to write CURL data to memory
size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* s) {
    size_t newLength = size * nmemb;
    try {
        s->append((char*)contents, newLength);
        return newLength;
    } catch(std::bad_alloc& e) {
        return 0;
    }
}

// Function to download map tile from LINZ
cv::Mat downloadMapTile(double north, double south, double east, double west, int width, int height) {
    CURL* curl;
    CURLcode res;
    std::string readBuffer;
    
    // Calculate the center of the bounding box
    double centerLat = (north + south) / 2.0;
    double centerLon = (east + west) / 2.0;
    
    // Dynamically calculate zoom level based on region width and output image width.
    // At zoom level z, one tile covers a normalized width of 1/(2^z) (with nominal tile size 256px).
    // Using our conversion, the relation is:
    //   2^z ≈ (width * 360) / (256 * (east - west))
    // Thus, z = log2((width * 360) / (256 * (east - west)))
    double regionWidthDeg = east - west; // in degrees
    double computed = (width * 360.0) / (256.0 * regionWidthDeg);
    int calculatedZoom = static_cast<int>(std::log2(computed));
    // Clamp zoom to a valid range (e.g., 0 to 19)
    if (calculatedZoom < 0) calculatedZoom = 0;
    if (calculatedZoom > 19) calculatedZoom = 19;
    int zoom = calculatedZoom;
    
    // Use LINZ basemap API with aerial imagery
    // Format: https://basemaps.linz.govt.nz/v1/tiles/{tileset_name}/{crs}/{z}/{x}/{y}.{format}?api={api_key}
    
    // Convert lat/lon to tile coordinates (Web Mercator)
    int x = static_cast<int>((centerLon + 180.0) / 360.0 * (1 << zoom));
    int y = static_cast<int>((1.0 - log(tan(centerLat * M_PI / 180.0) + 1.0 / cos(centerLat * M_PI / 180.0)) / M_PI) / 2.0 * (1 << zoom));
    y = (1 << zoom) - 1 - y; // Flip y coordinate for LINZ tile service to correct offset
    
    // Use aerial imagery with EPSG:3857 projection (Web Mercator)
    std::string url = "https://basemaps.linz.govt.nz/v1/tiles/aerial/EPSG:3857/" + 
                      std::to_string(zoom) + "/" + 
                      std::to_string(x) + "/" + 
                      std::to_string(y) + 
                      ".png?api=d01jrm3t2gzdycm5j8rh03e69fw";
    
    std::cout << "Downloading map from LINZ: " << url << std::endl;
    std::cout << "Note: Using LINZ demo API key. For production use, get your own key." << std::endl;
    
    curl = curl_easy_init();
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        
        // Set up headers to look like a browser request
        struct curl_slist *headers = NULL;
        headers = curl_slist_append(headers, "User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36");
        headers = curl_slist_append(headers, "Accept: image/jpeg,image/png,*/*");
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
        
        res = curl_easy_perform(curl);
        
        if(res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
            curl_easy_cleanup(curl);
            curl_slist_free_all(headers);
            
            // Try fallback to OpenStreetMap if LINZ fails
            return downloadOpenStreetMapTile(north, south, east, west, width, height);
        }
        
        curl_easy_cleanup(curl);
        curl_slist_free_all(headers);
        
        // Convert response to OpenCV image
        if (!readBuffer.empty()) {
            std::vector<char> data(readBuffer.begin(), readBuffer.end());
            cv::Mat img = cv::imdecode(data, cv::IMREAD_COLOR);
            
            if (!img.empty()) {
                // Resize to requested dimensions
                cv::resize(img, img, cv::Size(width, height));
                return img;
            } else {
                std::cerr << "Failed to decode image data from LINZ" << std::endl;
                // Try fallback
                return downloadOpenStreetMapTile(north, south, east, west, width, height);
            }
        } else {
            std::cerr << "Empty response from LINZ server" << std::endl;
            // Try fallback
            return downloadOpenStreetMapTile(north, south, east, west, width, height);
        }
    } else {
        std::cerr << "Failed to initialize curl" << std::endl;
    }
    
    // Return empty image if download fails
    return cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
}

// Fallback function to download map tile from OpenStreetMap
cv::Mat downloadOpenStreetMapTile(double north, double south, double east, double west, int width, int height) {
    CURL* curl;
    CURLcode res;
    std::string readBuffer;
    
    std::cout << "Falling back to OpenStreetMap..." << std::endl;
    
    // Calculate the center of the bounding box
    double centerLat = (north + south) / 2.0;
    double centerLon = (east + west) / 2.0;
    
    // Calculate appropriate zoom level based on the bounding box size
    double latDiff = north - south;
    double lonDiff = east - west;
    int zoom = 12; // Default zoom level
    
    // Adjust zoom level based on the area size
    if (latDiff > 0.5 || lonDiff > 0.5) zoom = 10;
    else if (latDiff < 0.05 || lonDiff < 0.05) zoom = 15;
    
    // Convert lat/lon to tile coordinates (Web Mercator)
    int x = static_cast<int>((centerLon + 180.0) / 360.0 * (1 << zoom));
    int y = static_cast<int>((1.0 - log(tan(centerLat * M_PI / 180.0) + 1.0 / cos(centerLat * M_PI / 180.0)) / M_PI) / 2.0 * (1 << zoom));
    
    // Use OpenStreetMap tile server
    std::string url = "https://tile.openstreetmap.org/" + 
                      std::to_string(zoom) + "/" + 
                      std::to_string(x) + "/" + 
                      std::to_string(y) + 
                      ".png";
    
    std::cout << "Downloading map from OpenStreetMap: " << url << std::endl;
    
    curl = curl_easy_init();
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        
        // Set up headers to look like a browser request
        struct curl_slist *headers = NULL;
        headers = curl_slist_append(headers, "User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36");
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
        
        res = curl_easy_perform(curl);
        
        if(res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
            curl_easy_cleanup(curl);
            curl_slist_free_all(headers);
            return cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
        }
        
        curl_easy_cleanup(curl);
        curl_slist_free_all(headers);
        
        // Convert response to OpenCV image
        if (!readBuffer.empty()) {
            std::vector<char> data(readBuffer.begin(), readBuffer.end());
            cv::Mat img = cv::imdecode(data, cv::IMREAD_COLOR);
            
            if (!img.empty()) {
                // Resize to requested dimensions
                cv::resize(img, img, cv::Size(width, height));
                return img;
            } else {
                std::cerr << "Failed to decode image data from OpenStreetMap" << std::endl;
            }
        } else {
            std::cerr << "Empty response from OpenStreetMap server" << std::endl;
        }
    }
    
    // If all attempts fail, create a blank image with grid lines
    cv::Mat blankImage(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
    
    // Draw grid lines
    for (int i = 0; i < width; i += 100) {
        cv::line(blankImage, cv::Point(i, 0), cv::Point(i, height), cv::Scalar(200, 200, 200), 1);
    }
    for (int i = 0; i < height; i += 100) {
        cv::line(blankImage, cv::Point(0, i), cv::Point(width, i), cv::Scalar(200, 200, 200), 1);
    }
    
    // Draw coordinate frame
    cv::putText(blankImage, "N: " + std::to_string(north), cv::Point(10, 20), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    cv::putText(blankImage, "S: " + std::to_string(south), cv::Point(10, height - 10), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    cv::putText(blankImage, "W: " + std::to_string(west), cv::Point(10, height / 2), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    cv::putText(blankImage, "E: " + std::to_string(east), cv::Point(width - 100, height / 2), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    
    return blankImage;
}

// Function to parse GPX file
std::vector<GpxPoint> parseGpxFile(const std::string& filename) {
    std::vector<GpxPoint> points;
    xmlDocPtr doc;
    xmlNodePtr cur;
    
    // Parse the GPX file
    doc = xmlParseFile(filename.c_str());
    if (doc == NULL) {
        std::cerr << "Error: could not parse file " << filename << std::endl;
        return points;
    }
    
    // Get the root element
    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
        std::cerr << "Error: empty document" << std::endl;
        xmlFreeDoc(doc);
        return points;
    }
    
    // Check if root element is gpx
    if (xmlStrcmp(cur->name, (const xmlChar *)"gpx")) {
        std::cerr << "Error: document of the wrong type, root node != gpx" << std::endl;
        xmlFreeDoc(doc);
        return points;
    }
    
    // Find trk element
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"trk"))) {
            xmlNodePtr trkChild = cur->xmlChildrenNode;
            while (trkChild != NULL) {
                if ((!xmlStrcmp(trkChild->name, (const xmlChar *)"trkseg"))) {
                    xmlNodePtr ptNode = trkChild->xmlChildrenNode;
                    while (ptNode != NULL) {
                        if ((!xmlStrcmp(ptNode->name, (const xmlChar *)"trkpt"))) {
                            GpxPoint point;
                            
                            // Get latitude and longitude attributes
                            xmlChar* lat = xmlGetProp(ptNode, (const xmlChar *)"lat");
                            xmlChar* lon = xmlGetProp(ptNode, (const xmlChar *)"lon");
                            
                            if (lat && lon) {
                                point.lat = std::stod((char*)lat);
                                point.lon = std::stod((char*)lon);
                                
                                // Get elevation and time if available
                                xmlNodePtr trkptChild = ptNode->xmlChildrenNode;
                                while (trkptChild != NULL) {
                                    if ((!xmlStrcmp(trkptChild->name, (const xmlChar *)"ele"))) {
                                        xmlChar* ele = xmlNodeListGetString(doc, trkptChild->xmlChildrenNode, 1);
                                        if (ele) {
                                            point.ele = std::stod((char*)ele);
                                            xmlFree(ele);
                                        }
                                    } else if ((!xmlStrcmp(trkptChild->name, (const xmlChar *)"time"))) {
                                        xmlChar* time = xmlNodeListGetString(doc, trkptChild->xmlChildrenNode, 1);
                                        if (time) {
                                            point.time = (char*)time;
                                            xmlFree(time);
                                        }
                                    }
                                    trkptChild = trkptChild->next;
                                }
                                
                                points.push_back(point);
                            }
                            
                            if (lat) xmlFree(lat);
                            if (lon) xmlFree(lon);
                        }
                        ptNode = ptNode->next;
                    }
                }
                trkChild = trkChild->next;
            }
        }
        cur = cur->next;
    }
    
    xmlFreeDoc(doc);
    return points;
}

// Function to calculate the bounding box of GPX points
void calculateBoundingBox(const std::vector<GpxPoint>& points, 
                         double& north, double& south, double& east, double& west) {
    if (points.empty()) return;
    
    north = south = points[0].lat;
    east = west = points[0].lon;
    
    for (const auto& point : points) {
        north = std::max(north, point.lat);
        south = std::min(south, point.lat);
        east = std::max(east, point.lon);
        west = std::min(west, point.lon);
    }
}

// Function to calculate the diagonal distance of the bounding box
double calculateDiagonalDistance(double north, double south, double east, double west) {
    // Approximate conversion from degrees to meters at the equator
    const double DEG_TO_METERS = 111319.9;
    
    double latDiff = (north - south) * DEG_TO_METERS;
    double lonDiff = (east - west) * DEG_TO_METERS * cos((north + south) * M_PI / 360.0);
    
    return sqrt(latDiff * latDiff + lonDiff * lonDiff);
}

// Function to calculate the distance between two points
double calculateDistance(const GpxPoint& p1, const GpxPoint& p2) {
    // Approximate conversion from degrees to meters at the equator
    const double DEG_TO_METERS = 111319.9;
    
    double latDiff = (p1.lat - p2.lat) * DEG_TO_METERS;
    double lonDiff = (p1.lon - p2.lon) * DEG_TO_METERS * cos((p1.lat + p2.lat) * M_PI / 360.0);
    
    return sqrt(latDiff * latDiff + lonDiff * lonDiff);
}

// Function to simplify GPX points
std::vector<GpxPoint> simplifyGpxPoints(const std::vector<GpxPoint>& points, double thresholdPercent) {
    if (points.size() <= 2) return points;
    
    std::vector<GpxPoint> simplifiedPoints;
    
    // Calculate bounding box
    double north, south, east, west;
    calculateBoundingBox(points, north, south, east, west);
    
    // Calculate diagonal distance and threshold
    double diagonal = calculateDiagonalDistance(north, south, east, west);
    double threshold = diagonal * thresholdPercent / 100.0;
    
    // Always keep the first point
    simplifiedPoints.push_back(points[0]);
    
    size_t i = 1;
    while (i < points.size() - 1) {
        std::vector<GpxPoint> cluster;
        cluster.push_back(points[i]);
        
        size_t j = i + 1;
        while (j < points.size() - 1 && calculateDistance(points[i], points[j]) < threshold) {
            cluster.push_back(points[j]);
            j++;
        }
        
        // Calculate median point for the cluster
        if (!cluster.empty()) {
            std::vector<double> lats, lons;
            for (const auto& p : cluster) {
                lats.push_back(p.lat);
                lons.push_back(p.lon);
            }
            
            std::sort(lats.begin(), lats.end());
            std::sort(lons.begin(), lons.end());
            
            GpxPoint medianPoint;
            medianPoint.lat = lats[lats.size() / 2];
            medianPoint.lon = lons[lons.size() / 2];
            
            // Use elevation and time from the middle point if available
            size_t midIndex = cluster.size() / 2;
            medianPoint.ele = cluster[midIndex].ele;
            medianPoint.time = cluster[midIndex].time;
            
            simplifiedPoints.push_back(medianPoint);
        }
        
        i = j;
    }
    
    // Always keep the last point
    simplifiedPoints.push_back(points.back());
    
    return simplifiedPoints;
}

// Function to write GPX file
void writeGpxFile(const std::string& filename, const std::vector<GpxPoint>& points) {
    xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
    xmlNodePtr root_node = xmlNewNode(NULL, (const xmlChar*)"gpx");
    
    // Add GPX attributes
    xmlNewProp(root_node, (const xmlChar*)"version", (const xmlChar*)"1.1");
    xmlNewProp(root_node, (const xmlChar*)"creator", (const xmlChar*)"gpxnicer");
    xmlNewProp(root_node, (const xmlChar*)"xmlns", (const xmlChar*)"http://www.topografix.com/GPX/1/1");
    
    xmlDocSetRootElement(doc, root_node);
    
    // Create track
    xmlNodePtr trk_node = xmlNewChild(root_node, NULL, (const xmlChar*)"trk", NULL);
    xmlNodePtr trkseg_node = xmlNewChild(trk_node, NULL, (const xmlChar*)"trkseg", NULL);
    
    // Add track points
    for (const auto& point : points) {
        xmlNodePtr trkpt_node = xmlNewChild(trkseg_node, NULL, (const xmlChar*)"trkpt", NULL);
        
        // Add latitude and longitude attributes
        xmlNewProp(trkpt_node, (const xmlChar*)"lat", (const xmlChar*)std::to_string(point.lat).c_str());
        xmlNewProp(trkpt_node, (const xmlChar*)"lon", (const xmlChar*)std::to_string(point.lon).c_str());
        
        // Add elevation if available
        if (point.ele != 0.0) {
            xmlNewChild(trkpt_node, NULL, (const xmlChar*)"ele", (const xmlChar*)std::to_string(point.ele).c_str());
        }
        
        // Add time if available
        if (!point.time.empty()) {
            xmlNewChild(trkpt_node, NULL, (const xmlChar*)"time", (const xmlChar*)point.time.c_str());
        }
    }
    
    // Save the document
    xmlSaveFormatFileEnc(filename.c_str(), doc, "UTF-8", 1);
    xmlFreeDoc(doc);
}

// Function to create visualization
void createVisualization(const std::string& filename, 
                        const std::vector<GpxPoint>& originalPoints,
                        const std::vector<GpxPoint>& simplifiedPoints) {
    // Calculate bounding box
    double north, south, east, west;
    calculateBoundingBox(originalPoints, north, south, east, west);
    
    // Add some padding to the bounding box
    double latPadding = (north - south) * 0.1;
    double lonPadding = (east - west) * 0.1;
    north += latPadding;
    south -= latPadding;
    east += lonPadding;
    west -= lonPadding;

    // Calculate Web Mercator bounds from the padded lat/lon box
    double mercatorNW_x, mercatorNW_y, mercatorSE_x, mercatorSE_y;
    latLonToWebMercator(north, west, mercatorNW_x, mercatorNW_y);
    latLonToWebMercator(south, east, mercatorSE_x, mercatorSE_y);
    double xPadding = (mercatorSE_x - mercatorNW_x) * 0.05;
    double yPadding = (mercatorSE_y - mercatorNW_y) * 0.05;
    double westMerc = mercatorNW_x - xPadding;
    double eastMerc = mercatorSE_x + xPadding;
    double northMerc = mercatorNW_y - yPadding;
    double southMerc = mercatorSE_y + yPadding;

    // Create image
    int width = 1200;
    int height = 800;
    cv::Mat img = downloadMapTile(north, south, east, west, width, height);
    
    // Check if the downloaded image is valid, if not create a blank image
    if (img.empty()) {
        std::cerr << "Warning: Could not download map tile. Creating blank image instead." << std::endl;
        img = cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
    }
    
    // Calculate Web Mercator for all four corners
    double nw_x, nw_y, ne_x, ne_y, sw_x, sw_y, se_x, se_y;
    latLonToWebMercator(north, west, nw_x, nw_y);
    latLonToWebMercator(north, east, ne_x, ne_y);
    latLonToWebMercator(south, west, sw_x, sw_y);
    latLonToWebMercator(south, east, se_x, se_y);

    // Debug output
    std::cout << "Bounding box in Web Mercator:" << std::endl;
    std::cout << "North: " << northMerc << ", South: " << southMerc << std::endl;
    std::cout << "West: " << westMerc << ", East: " << eastMerc << std::endl;
    
    // Draw original track in blue
    for (size_t i = 1; i < originalPoints.size(); i++) {
        double x1, y1, x2, y2;
        latLonToWebMercator(originalPoints[i-1].lat, originalPoints[i-1].lon, x1, y1);
        latLonToWebMercator(originalPoints[i].lat, originalPoints[i].lon, x2, y2);
        
        cv::Point p1 = webMercatorToPixel(x1, y1, westMerc, eastMerc, northMerc, southMerc, width, height);
        cv::Point p2 = webMercatorToPixel(x2, y2, westMerc, eastMerc, northMerc, southMerc, width, height);
        
        // Debug output for first few points
        if (i < 5) {
            std::cout << "Point " << i << ":" << std::endl;
            std::cout << "  Lat/Lon: " << originalPoints[i].lat << ", " << originalPoints[i].lon << std::endl;
            std::cout << "  Web Mercator: " << x2 << ", " << y2 << std::endl;
            std::cout << "  Pixel: " << p2.x << ", " << p2.y << std::endl;
        }
        
        // Only draw if points are within image bounds
        if (p1.x >= 0 && p1.x < width && p1.y >= 0 && p1.y < height &&
            p2.x >= 0 && p2.x < width && p2.y >= 0 && p2.y < height) {
            cv::line(img, p1, p2, cv::Scalar(255, 0, 0), 2);
        }
    }
    
    // Draw simplified track in red
    for (size_t i = 1; i < simplifiedPoints.size(); i++) {
        double x1, y1, x2, y2;
        latLonToWebMercator(simplifiedPoints[i-1].lat, simplifiedPoints[i-1].lon, x1, y1);
        latLonToWebMercator(simplifiedPoints[i].lat, simplifiedPoints[i].lon, x2, y2);
        
        cv::Point p1 = webMercatorToPixel(x1, y1, westMerc, eastMerc, northMerc, southMerc, width, height);
        cv::Point p2 = webMercatorToPixel(x2, y2, westMerc, eastMerc, northMerc, southMerc, width, height);
        
        // Only draw if points are within image bounds
        if (p1.x >= 0 && p1.x < width && p1.y >= 0 && p1.y < height &&
            p2.x >= 0 && p2.x < width && p2.y >= 0 && p2.y < height) {
            cv::line(img, p1, p2, cv::Scalar(0, 0, 255), 2);
        }
    }
    
    // Add legend
    cv::rectangle(img, cv::Point(10, 10), cv::Point(200, 60), cv::Scalar(255, 255, 255), -1);
    cv::line(img, cv::Point(20, 25), cv::Point(50, 25), cv::Scalar(255, 0, 0), 2);
    cv::putText(img, "Original", cv::Point(60, 30), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    cv::line(img, cv::Point(20, 45), cv::Point(50, 45), cv::Scalar(0, 0, 255), 2);
    cv::putText(img, "Simplified", cv::Point(60, 50), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    
    // Add copyright notice
    cv::rectangle(img, cv::Point(10, height - 30), cv::Point(350, height - 10), cv::Scalar(255, 255, 255, 128), -1);
    cv::putText(img, "Map data: LINZ / OpenStreetMap contributors", cv::Point(15, height - 15), 
                cv::FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(0, 0, 0), 1);
    
    // Save image
    try {
        bool success = cv::imwrite(filename, img);
        if (!success) {
            std::cerr << "Error: Failed to write image to " << filename << std::endl;
        }
    } catch (const cv::Exception& ex) {
        std::cerr << "Exception while saving image: " << ex.what() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <gpx_file> <threshold_percent>" << std::endl;
        return 1;
    }
    
    std::string inputFile = argv[1];
    double thresholdPercent;
    
    try {
        thresholdPercent = std::stod(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << "Error: threshold_percent must be a valid number" << std::endl;
        return 1;
    }
    
    // Check if threshold is valid
    if (thresholdPercent <= 0 || thresholdPercent >= 100) {
        std::cerr << "Error: threshold_percent must be between 0 and 100" << std::endl;
        return 1;
    }
    
    // Check if input file exists
    std::ifstream fileCheck(inputFile);
    if (!fileCheck.good()) {
        std::cerr << "Error: cannot open file " << inputFile << std::endl;
        return 1;
    }
    fileCheck.close();
    
    // Parse GPX file
    std::vector<GpxPoint> originalPoints = parseGpxFile(inputFile);
    if (originalPoints.empty()) {
        std::cerr << "Error: no points found in GPX file" << std::endl;
        return 1;
    }
    
    // Simplify points
    std::vector<GpxPoint> simplifiedPoints = simplifyGpxPoints(originalPoints, thresholdPercent);
    
    // Generate output filenames
    std::string outputGpxFile = inputFile.substr(0, inputFile.find_last_of('.')) + "_simplified.gpx";
    std::string outputJpgFile = inputFile.substr(0, inputFile.find_last_of('.')) + "_simplified.jpg";
    
    // Write simplified GPX file
    try {
        writeGpxFile(outputGpxFile, simplifiedPoints);
        std::cout << "Successfully wrote GPX file: " << outputGpxFile << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error writing GPX file: " << e.what() << std::endl;
        return 1;
    }
    
    // Create visualization
    try {
        createVisualization(outputJpgFile, originalPoints, simplifiedPoints);
        std::cout << "Successfully created visualization: " << outputJpgFile << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error creating visualization: " << e.what() << std::endl;
        // Continue execution even if visualization fails
    }
    
    // Print summary
    std::cout << "Original points: " << originalPoints.size() << std::endl;
    std::cout << "Simplified points: " << simplifiedPoints.size() << std::endl;
    std::cout << "Reduction: " << (1.0 - static_cast<double>(simplifiedPoints.size()) / originalPoints.size()) * 100.0 << "%" << std::endl;
    std::cout << "Output GPX file: " << outputGpxFile << std::endl;
    std::cout << "Output JPG file: " << outputJpgFile << std::endl;
    
    return 0;
}
