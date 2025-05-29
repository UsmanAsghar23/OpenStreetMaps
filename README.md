# 🗺️ OpenStreetMaps

This C++ project allows users to navigate the University of Illinois at Chicago (UIC) campus using real OpenStreetMap (OSM) data. It models walking paths and buildings, and uses Dijkstra's algorithm to find the shortest walking route between two people’s buildings and a mutually close destination.

## 📌 Features

- Load and parse OSM files using TinyXML2.
- Read and process:
  - Map nodes (lat/lon coordinates)
  - Footways (walkable paths)
  - Campus buildings (name, abbreviation, location)
- Build a graph from footways.
- Use Dijkstra’s algorithm to compute shortest paths.
- Interactively prompt two users for their buildings and compute:
  - Their nearest walkable nodes
  - Shortest walking path to a shared destination (midpoint)
  - Distance and node-by-node path for each person

## 📂 Project Structure

| File               | Description |
|--------------------|-------------|
| `application.cpp`  | Main program logic and user interaction |
| `graph.h`          | Templated graph implementation using an adjacency list |
| `dist.cpp/h`       | Calculates distance and midpoint between two GPS coordinates |
| `osm.cpp/h`        | Parses `.osm` XML files and extracts relevant map data |
| `tinyxml2.cpp/h`   | XML parser library used to read OSM files |
| `graph.txt`        | Sample graph input for testing |
| `graph.pdf`        | Visual representation of the sample graph |
| `map.osm` / `uic.osm` | XML files containing OpenStreetMap data |
| `testing.cpp`      | Standalone tester for `graph.h` using `graph.txt` |
| `makefile`         | Build file for compiling the project |

## 🛠️ Building the Project

### ✅ Prerequisites

- C++17 or later
- A compiler like `g++`
- Make (if using the provided `makefile`)

### 🔧 Build with `makefile`

```bash
make
