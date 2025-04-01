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

// Convert WGS84 (lat/lon) to Web Mercator coordinates
void latLonToWebMercator(double lat, double lon, double& x, double& y) {
    // Convert to radians
    lat = lat * M_PI / 180.0;
    lon = lon * M_PI / 180.0;
    
    // Web Mercator projection (EPSG:3857)
    x = (lon + M_PI) / (2.0 * M_PI);
    y = (1.0 - log(tan(lat) + 1.0 / cos(lat)) / M_PI) / 2.0;
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

// Function to download satellite imagery from LINZ Basemaps API
cv::Mat downloadLinzSatelliteImage(double north, double south, double east, double west, int width, int height) {
    CURL* curl;
    CURLcode res;
    
    // Create a blank canvas to build the final image
    cv::Mat finalImage(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
    
    // Calculate proper zoom level for the area
    // This formula gives better results when we need to see details of a track
    double latDiff = north - south;
    double lonDiff = east - west;
    double maxDiff = std::max(latDiff, lonDiff);
    
    // Adjust zoom based on the geographic area size
    // The smaller the area, the higher the zoom level
    int zoom = 15; // Start with a reasonable default for GPX tracks
    
    if (maxDiff > 1.0) zoom = 9;      // Very large area
    else if (maxDiff > 0.5) zoom = 10; // Large area
    else if (maxDiff > 0.1) zoom = 12; // Medium area
    else if (maxDiff > 0.05) zoom = 13; // Medium-small area
    else if (maxDiff > 0.01) zoom = 14; // Small area
    else if (maxDiff > 0.005) zoom = 15; // Very small area
    else zoom = 16;                    // Extremely small area
    
    // Ensure zoom is within valid range (0-22 for Web Mercator)
    if (zoom < 0) zoom = 0;
    if (zoom > 20) zoom = 20; // Cap at 20 for safety
    
    std::cout << "Using zoom level: " << zoom << " for area width: " << maxDiff << " degrees" << std::endl;
    
    // LINZ Basemap API key (standard access)
    const std::string api_key = "d01jrm3t2gzdycm5j8rh03e69fw";  // Demo key from LINZ docs
    
    // Tile to lat/lon conversion helpers for this zoom level
    auto tileToLon = [zoom](int x) -> double {
        return x / std::pow(2.0, zoom) * 360.0 - 180.0;
    };
    
    auto tileToLat = [zoom](int y) -> double {
        double n = M_PI - 2.0 * M_PI * y / std::pow(2.0, zoom);
        return 180.0 / M_PI * std::atan(0.5 * (std::exp(n) - std::exp(-n)));
    };
    
    // Function to convert lat/lon to tile numbers
    auto latLonToTile = [zoom](double lat, double lon, int& x, int& y) {
        // Convert to radians
        double latRad = lat * M_PI / 180.0;
        
        // Calculate tile coordinates
        x = static_cast<int>((lon + 180.0) / 360.0 * (1 << zoom));
        y = static_cast<int>((1.0 - log(tan(latRad) + 1.0 / cos(latRad)) / M_PI) / 2.0 * (1 << zoom));
    };
    
    // Calculate tile range that covers our area
    int minX, minY, maxX, maxY;
    latLonToTile(north, west, minX, minY);
    latLonToTile(south, east, maxX, maxY);
    
    // Ensure we have the correct order
    if (minX > maxX) std::swap(minX, maxX);
    if (minY > maxY) std::swap(minY, maxY);
    
    // Limit the number of tiles to prevent excessive requests
    const int MAX_TILES = 9; // 3x3 grid max
    int tilesX = maxX - minX + 1;
    int tilesY = maxY - minY + 1;
    
    if (tilesX * tilesY > MAX_TILES) {
        // If too many tiles, adjust to the center of the area
        int centerX = (minX + maxX) / 2;
        int centerY = (minY + maxY) / 2;
        
        // This gives us a 3x3 grid of tiles centered on the area of interest
        minX = centerX - 1;
        maxX = centerX + 1;
        minY = centerY - 1;
        maxY = centerY + 1;
        
        tilesX = maxX - minX + 1;
        tilesY = maxY - minY + 1;
    }
    
    // Calculate the exact lat/lon bounds of this tile grid
    double tileGridWest = tileToLon(minX);
    double tileGridEast = tileToLon(maxX + 1);  // +1 because tiles cover [minX, maxX+1)
    double tileGridNorth = tileToLat(minY);
    double tileGridSouth = tileToLat(maxY + 1); // +1 because tiles cover [minY, maxY+1)
    
    std::cout << "Downloading " << tilesX << "x" << tilesY << " tiles from LINZ" << std::endl;
    std::cout << "Tile grid covers: North=" << tileGridNorth << ", South=" << tileGridSouth
              << ", East=" << tileGridEast << ", West=" << tileGridWest << std::endl;
    std::cout << "GPX area covers: North=" << north << ", South=" << south
              << ", East=" << east << ", West=" << west << std::endl;
    
    // Download and stitch tiles
    std::vector<cv::Mat> tiles;
    std::vector<std::pair<int, int>> tilePositions;
    const int tileSize = 256; // Standard tile size
    
    curl = curl_easy_init();
    if (!curl) {
        std::cerr << "Failed to initialize curl" << std::endl;
        return finalImage;
    }
    
    // Set up headers to look like a browser request
    struct curl_slist *headers = NULL;
    headers = curl_slist_append(headers, "User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    
    // Try to download a single tile first to test connectivity
    std::string testBuffer;
    int testX = (minX + maxX) / 2;
    int testY = (minY + maxY) / 2;
    
    // Build the LINZ Basemap URL with PNG format
    std::string testUrl = "https://basemaps.linz.govt.nz/v1/tiles/aerial/EPSG:3857/" + 
                          std::to_string(zoom) + "/" + 
                          std::to_string(testX) + "/" + 
                          std::to_string(testY) + 
                          ".png?api=" + api_key;  // Use PNG format which is more widely supported
    
    std::cout << "Testing connection with tile URL: " << testUrl << std::endl;
    
    curl_easy_setopt(curl, CURLOPT_URL, testUrl.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &testBuffer);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L); // 10 second timeout
    
    res = curl_easy_perform(curl);
    
    if (res != CURLE_OK) {
        std::cerr << "Test tile download failed: " << curl_easy_strerror(res) << std::endl;
        std::cerr << "Check your internet connection or the LINZ API status." << std::endl;
        curl_slist_free_all(headers);
        curl_easy_cleanup(curl);
        return finalImage;
    }
    
    // Try to decode the test image
    bool canUseLinz = false;
    if (!testBuffer.empty()) {
        std::vector<char> testData(testBuffer.begin(), testBuffer.end());
        cv::Mat testImg = cv::imdecode(testData, cv::IMREAD_COLOR);
        canUseLinz = !testImg.empty();
        
        if (canUseLinz) {
            std::cout << "Successfully connected to LINZ API. Downloading tiles..." << std::endl;
        } else {
            std::cerr << "Failed to decode test tile. Falling back to blank background." << std::endl;
            curl_slist_free_all(headers);
            curl_easy_cleanup(curl);
            return finalImage;
        }
    } else {
        std::cerr << "Empty response from LINZ server for test tile." << std::endl;
        curl_slist_free_all(headers);
        curl_easy_cleanup(curl);
        return finalImage;
    }
    
    // Download each tile
    bool anyTileSuccess = false;
    for (int y = minY; y <= maxY; y++) {
        for (int x = minX; x <= maxX; x++) {
            std::string readBuffer;
            
            // Build the LINZ Basemap URL with PNG format
            std::string url = "https://basemaps.linz.govt.nz/v1/tiles/aerial/EPSG:3857/" + 
                             std::to_string(zoom) + "/" + 
                             std::to_string(x) + "/" + 
                             std::to_string(y) + 
                             ".png?api=" + api_key;  // Use PNG format
            
            if (x == minX && y == minY) {
                std::cout << "Downloading tiles with URL pattern: " << url << std::endl;
                std::cout << "Note: Using LINZ demo API key. For production use, get your own key." << std::endl;
            }
            
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
            
            res = curl_easy_perform(curl);
            
            if (res != CURLE_OK) {
                std::cerr << "Failed to download tile at x=" << x << ", y=" << y << ": " << curl_easy_strerror(res) << std::endl;
                continue;
            }
            
            // Convert response to OpenCV image
            if (!readBuffer.empty()) {
                std::vector<char> data(readBuffer.begin(), readBuffer.end());
                cv::Mat tileImg = cv::imdecode(data, cv::IMREAD_COLOR);
                
                if (!tileImg.empty()) {
                    tiles.push_back(tileImg);
                    tilePositions.push_back(std::make_pair(x - minX, y - minY));
                    anyTileSuccess = true;
                } else {
                    std::cerr << "Failed to decode tile image at x=" << x << ", y=" << y << std::endl;
                }
            }
        }
    }
    
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
    
    if (!anyTileSuccess) {
        std::cerr << "Failed to download any tiles from LINZ" << std::endl;
        return finalImage;
    }
    
    // Create a larger canvas to hold all tiles
    int fullWidth = tilesX * tileSize;
    int fullHeight = tilesY * tileSize;
    cv::Mat stitchedImage(fullHeight, fullWidth, CV_8UC3, cv::Scalar(255, 255, 255));
    
    // Place tiles in the correct positions
    for (size_t i = 0; i < tiles.size(); i++) {
        int tileX = tilePositions[i].first;
        int tileY = tilePositions[i].second;
        
        cv::Rect roi(tileX * tileSize, tileY * tileSize, tileSize, tileSize);
        if (roi.x + roi.width <= stitchedImage.cols && roi.y + roi.height <= stitchedImage.rows) {
            tiles[i].copyTo(stitchedImage(roi));
        }
    }
    
    // Now we need to crop the stitched image to match exactly our GPX bounds
    // Find the pixel coordinates in the stitched image that correspond to our GPX bounds
    auto lonToPixelX = [tileGridWest, tileGridEast, fullWidth](double lon) -> int {
        return static_cast<int>((lon - tileGridWest) / (tileGridEast - tileGridWest) * fullWidth);
    };
    
    auto latToPixelY = [tileGridNorth, tileGridSouth, fullHeight](double lat) -> int {
        return static_cast<int>((tileGridNorth - lat) / (tileGridNorth - tileGridSouth) * fullHeight);
    };
    
    // Calculate crop region
    int cropX = lonToPixelX(west);
    int cropY = latToPixelY(north);
    int cropWidth = lonToPixelX(east) - cropX;
    int cropHeight = latToPixelY(south) - cropY;
    
    // Ensure crop region is valid
    cropX = std::max(0, std::min(fullWidth - 1, cropX));
    cropY = std::max(0, std::min(fullHeight - 1, cropY));
    cropWidth = std::max(1, std::min(fullWidth - cropX, cropWidth));
    cropHeight = std::max(1, std::min(fullHeight - cropY, cropHeight));
    
    // Crop the image
    cv::Mat croppedImage;
    cv::Rect cropROI(cropX, cropY, cropWidth, cropHeight);
    croppedImage = stitchedImage(cropROI);
    
    // Resize to desired dimensions
    cv::resize(croppedImage, finalImage, cv::Size(width, height));
    
    // Apply a slight blur to make the satellite imagery less distracting
    cv::GaussianBlur(finalImage, finalImage, cv::Size(3, 3), 0);
    
    return finalImage;
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

// Function to create visualization with LINZ satellite imagery background
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

    // Create image
    int width = 1200;
    int height = 800;
    
    // Download satellite imagery from LINZ as background
    cv::Mat img = downloadLinzSatelliteImage(north, south, east, west, width, height);
    
    // Check if the downloaded image is valid, if not create a blank image with grid
    if (img.empty()) {
        std::cerr << "Warning: Could not download LINZ satellite imagery. Creating blank image instead." << std::endl;
        img = cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
        
        // Draw grid lines on blank background
        for (int i = 0; i < width; i += 100) {
            cv::line(img, cv::Point(i, 0), cv::Point(i, height), cv::Scalar(220, 220, 220), 1);
        }
        for (int i = 0; i < height; i += 100) {
            cv::line(img, cv::Point(0, i), cv::Point(width, i), cv::Scalar(220, 220, 220), 1);
        }
        
        // Draw coordinate labels
        cv::putText(img, "N: " + std::to_string(north), cv::Point(10, 20), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
        cv::putText(img, "S: " + std::to_string(south), cv::Point(10, height - 10), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
        cv::putText(img, "W: " + std::to_string(west), cv::Point(10, height / 2), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
        cv::putText(img, "E: " + std::to_string(east), cv::Point(width - 100, height / 2), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    } else {
        // Add subtle grid lines on satellite image for better orientation
        for (int i = 0; i < width; i += 200) {
            cv::line(img, cv::Point(i, 0), cv::Point(i, height), cv::Scalar(200, 200, 200, 128), 1);
        }
        for (int i = 0; i < height; i += 200) {
            cv::line(img, cv::Point(0, i), cv::Point(width, i), cv::Scalar(200, 200, 200, 128), 1);
        }
    }
    
    // Helper function to convert lat/lon to pixel coordinates
    auto latLonToPixel = [&](double lat, double lon) -> cv::Point {
        int x = static_cast<int>((lon - west) / (east - west) * width);
        int y = static_cast<int>((north - lat) / (north - south) * height);
        return cv::Point(x, y);
    };
    
    // Draw original track in blue
    for (size_t i = 1; i < originalPoints.size(); i++) {
        cv::Point p1 = latLonToPixel(originalPoints[i-1].lat, originalPoints[i-1].lon);
        cv::Point p2 = latLonToPixel(originalPoints[i].lat, originalPoints[i].lon);
        
        // Only draw if points are within image bounds
        if (p1.x >= 0 && p1.x < width && p1.y >= 0 && p1.y < height &&
            p2.x >= 0 && p2.x < width && p2.y >= 0 && p2.y < height) {
            cv::line(img, p1, p2, cv::Scalar(255, 0, 0), 2);
        }
    }
    
    // Draw simplified track in red
    for (size_t i = 1; i < simplifiedPoints.size(); i++) {
        cv::Point p1 = latLonToPixel(simplifiedPoints[i-1].lat, simplifiedPoints[i-1].lon);
        cv::Point p2 = latLonToPixel(simplifiedPoints[i].lat, simplifiedPoints[i].lon);
        
        // Only draw if points are within image bounds
        if (p1.x >= 0 && p1.x < width && p1.y >= 0 && p1.y < height &&
            p2.x >= 0 && p2.x < width && p2.y >= 0 && p2.y < height) {
            cv::line(img, p1, p2, cv::Scalar(0, 0, 255), 2);
        }
    }
    
    // Add legend with white background for visibility
    cv::rectangle(img, cv::Point(10, 10), cv::Point(200, 60), cv::Scalar(255, 255, 255, 200), -1);
    cv::line(img, cv::Point(20, 25), cv::Point(50, 25), cv::Scalar(255, 0, 0), 2);
    cv::putText(img, "Original", cv::Point(60, 30), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    cv::line(img, cv::Point(20, 45), cv::Point(50, 45), cv::Scalar(0, 0, 255), 2);
    cv::putText(img, "Simplified", cv::Point(60, 50), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    
    // Add attribution for LINZ data
    cv::rectangle(img, cv::Point(10, height - 30), cv::Point(350, height - 10), cv::Scalar(255, 255, 255, 200), -1);
    cv::putText(img, "Map data: LINZ Aerial Imagery", cv::Point(15, height - 15), 
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
