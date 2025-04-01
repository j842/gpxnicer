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

// Function to download map tile from OpenStreetMap
cv::Mat downloadMapTile(double north, double south, double east, double west, int width, int height) {
    CURL* curl;
    CURLcode res;
    std::string readBuffer;
    
    // Create URL for static map API (using OpenStreetMap Static Map API)
    std::string url = "https://maps.googleapis.com/maps/api/staticmap?center=" + 
                      std::to_string((north + south) / 2.0) + "," + 
                      std::to_string((east + west) / 2.0) + 
                      "&zoom=12&size=" + std::to_string(width) + "x" + std::to_string(height) + 
                      "&maptype=satellite&key=YOUR_API_KEY"; // Replace with your API key
    
    curl = curl_easy_init();
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        
        // Convert response to OpenCV image
        std::vector<char> data(readBuffer.begin(), readBuffer.end());
        cv::Mat img = cv::imdecode(data, cv::IMREAD_COLOR);
        return img;
    }
    
    // Return empty image if download fails
    return cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
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
    
    // Create image
    int width = 1200;
    int height = 800;
    cv::Mat img = downloadMapTile(north, south, east, west, width, height);
    
    // Function to convert lat/lon to pixel coordinates
    auto latLonToPixel = [&](double lat, double lon) -> cv::Point {
        int x = static_cast<int>((lon - west) / (east - west) * width);
        int y = static_cast<int>((north - lat) / (north - south) * height);
        return cv::Point(x, y);
    };
    
    // Draw original track in blue
    for (size_t i = 1; i < originalPoints.size(); i++) {
        cv::Point p1 = latLonToPixel(originalPoints[i-1].lat, originalPoints[i-1].lon);
        cv::Point p2 = latLonToPixel(originalPoints[i].lat, originalPoints[i].lon);
        cv::line(img, p1, p2, cv::Scalar(255, 0, 0), 2);
    }
    
    // Draw simplified track in red
    for (size_t i = 1; i < simplifiedPoints.size(); i++) {
        cv::Point p1 = latLonToPixel(simplifiedPoints[i-1].lat, simplifiedPoints[i-1].lon);
        cv::Point p2 = latLonToPixel(simplifiedPoints[i].lat, simplifiedPoints[i].lon);
        cv::line(img, p1, p2, cv::Scalar(0, 0, 255), 2);
    }
    
    // Add legend
    cv::rectangle(img, cv::Point(10, 10), cv::Point(200, 60), cv::Scalar(255, 255, 255), -1);
    cv::line(img, cv::Point(20, 25), cv::Point(50, 25), cv::Scalar(255, 0, 0), 2);
    cv::putText(img, "Original", cv::Point(60, 30), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    cv::line(img, cv::Point(20, 45), cv::Point(50, 45), cv::Scalar(0, 0, 255), 2);
    cv::putText(img, "Simplified", cv::Point(60, 50), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    
    // Save image
    cv::imwrite(filename, img);
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <gpx_file> <threshold_percent>" << std::endl;
        return 1;
    }
    
    std::string inputFile = argv[1];
    double thresholdPercent = std::stod(argv[2]);
    
    // Check if threshold is valid
    if (thresholdPercent <= 0 || thresholdPercent >= 100) {
        std::cerr << "Error: threshold_percent must be between 0 and 100" << std::endl;
        return 1;
    }
    
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
    writeGpxFile(outputGpxFile, simplifiedPoints);
    
    // Create visualization
    createVisualization(outputJpgFile, originalPoints, simplifiedPoints);
    
    // Print summary
    std::cout << "Original points: " << originalPoints.size() << std::endl;
    std::cout << "Simplified points: " << simplifiedPoints.size() << std::endl;
    std::cout << "Reduction: " << (1.0 - static_cast<double>(simplifiedPoints.size()) / originalPoints.size()) * 100.0 << "%" << std::endl;
    std::cout << "Output GPX file: " << outputGpxFile << std::endl;
    std::cout << "Output JPG file: " << outputJpgFile << std::endl;
    
    return 0;
}
