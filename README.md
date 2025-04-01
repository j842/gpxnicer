# gpxnicer

A C++ tool to simplify GPX files by reducing the number of points while preserving the overall path shape.

## Features

- Simplifies GPX files based on a configurable distance threshold
- Preserves the first and last points of the track
- Uses median latitude and longitude for simplified points
- Generates a visualization of the original and simplified tracks on a map
- Outputs both a simplified GPX file and a JPG visualization

## Usage

```
./gpxnicer input.gpx threshold_percent
```

Where:
- `input.gpx` is the path to the GPX file to simplify
- `threshold_percent` is the percentage of the bounding box diagonal to use as the simplification threshold

## Output

- A simplified GPX file named `input_simplified.gpx`
- A JPG visualization named `input_simplified.jpg` showing both the original (blue) and simplified (red) tracks on a map

## Building

```
mkdir build
cd build
cmake ..
make
```

## Dependencies

- C++17 or later
- libcurl
- libxml2
- OpenCV

## Map Data

This application uses map data from:
- LINZ (Land Information New Zealand) aerial imagery using their API
- If LINZ API access fails, a blank grid with coordinates is generated

Map data © LINZ, available under Creative Commons 4.0 International license.
