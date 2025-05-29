// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <stack>
#include <queue>
#include <algorithm>
#include <limits>

#include "tinyxml2.h"
#include "dist.h"
#include "graph.h"
#include "osm.h"

using namespace std;
using namespace tinyxml2;
const double INFINITY_VAL = numeric_limits<double>::max();

struct prioritize {
  //Overload operator used to compare pairs
  bool operator()(const pair<long long, double>& p1, const pair<long long, double>& p2) const {
    return p1.second > p2.second; //Compares the second element of the pair in decending order
  }
};

//Function to search for building based on name or abbreviation 
BuildingInfo searchForBuilding(const string& search, const vector<BuildingInfo>& Buildings) {
  for (const auto& building : Buildings) { //Iterate though each building in the vector that stores the buildings
    if (search == building.Abbrev || building.Fullname.find(search) != string::npos) { //Check if the search matches the name or abbreviaton 
      return building;
    }
  }
  return BuildingInfo(); //If no match is found return an empty building
}

//Find the closest building to a given midpoint based on distance
BuildingInfo findClosetBuilding(const Coordinates& midpoint, const vector<BuildingInfo>& buildings) {
  BuildingInfo nearestBuilding; //Varaible that stores information about the nearing buidling
  double minDist = INFINITY_VAL; //Initialize the minimum distace to high value
  for (const auto& building : buildings) { 
    //Calculate the distance between the midpoint and the current building
    double distance = distBetween2Points(midpoint.Lat, midpoint.Lon, building.Coords.Lat, building.Coords.Lon);
    if (distance < minDist) { //Check the calculated distance is less the than the current minimum distance
        minDist = distance;
        nearestBuilding = building; //Update the nearestBuilding information to the current building
    }
  }
  return nearestBuilding;
}

//This function finds the shortest paths from a starting node to all other nodes in the graph using Dijkstra algorithm.
void FindShortestPaths(graph<long long, double>& G, long long& startingNode,
                       map<long long, long long>& previousNodes,
                       map<long long, double>& nodeDistances) {
  //Priority queue to store nodes with their distances
  priority_queue<pair<long long, double>, vector<pair<long long, double>>, prioritize> pq;
  set<long long> visited; //Se to keep track of the nodes that have been visited
  for (auto& currentNode : G.getVertices()) { //Initialize distances and previous nodes for each node in the graph
    previousNodes[currentNode] = -1;
    nodeDistances[currentNode] = INFINITY_VAL;
    pq.push(make_pair(currentNode, INFINITY_VAL));
  }
  nodeDistances[startingNode] = 0; //Set distance of starting node to 0
  pq.push(make_pair(startingNode, 0));
  while (!pq.empty()) {
    pair<long long, double> currentNode = pq.top(); //Get the node with the smallest distance
    pq.pop();
    if (currentNode.second == INFINITY_VAL) { //Condtion breaks if remaining nodes are unreachable
      break;
    } 
    else if (visited.count(currentNode.first) >= 1) { //Is skipped if the node has been visited
      continue;
    } 
    else { //Marks the node as visited
      visited.insert(currentNode.first);
    }
    for (auto& adjacentNode : G.neighbors(currentNode.first)) { //Find the neighbors if the current Node
      double weight = 0;
      G.getWeight(currentNode.first, adjacentNode, weight);
      double otherPath = nodeDistances[currentNode.first] + weight; //Calculates the total distance from the current Node to the neighbor node
      if (otherPath < nodeDistances[adjacentNode]) { //Sees if a shorter path is found and updates the previous node
        nodeDistances[adjacentNode] = otherPath;
        pq.push(make_pair(adjacentNode, otherPath));
        previousNodes[adjacentNode] = currentNode.first;
      }
    }
  }
}

//Calculates and outputs the stats for navigating the UIC street map
void application(
  map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
  vector<BuildingInfo>& Buildings, graph<long long, double> G) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);
  //Main loop that handles everything
  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    //Search for building information for person 1 and person 2
    BuildingInfo abbrev1 = searchForBuilding(person1Building, Buildings);
    BuildingInfo abbrev2 = searchForBuilding(person2Building, Buildings);
    //Loops until valid buildings are entered for both people
    while (abbrev1.Fullname.empty() || abbrev2.Fullname.empty()) {
      if (abbrev1.Abbrev.empty()) {
        cout << "Person 1's building not found" << endl;
      }
      if (abbrev2.Abbrev.empty()) {
        cout << "Person 2's building not found" << endl;
      }
      cout << endl;
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      if (person1Building == "#") {
        return;
      }
      cout << "Enter person 2's building (partial name or abbreviation)> ";
      getline(cin, person2Building);
      abbrev1 = searchForBuilding(person1Building, Buildings);
      abbrev2 = searchForBuilding(person2Building, Buildings);
    }
    //Inintialize maps to store shortest paths and distances
    map<long long, long long> predOne;
    map<long long, long long> predTwo;
    map<long long, double> distanceOne;
    map<long long, double> distanceTwo;
    
    cout << endl;
    cout << "Person 1's point: " << endl;
    cout << " " << abbrev1.Fullname << endl;
    cout << " (" << abbrev1.Coords.Lat << ", " << abbrev1.Coords.Lon << ")" << endl;
    cout << "Person 2's point: " << endl;
    cout << " " << abbrev2.Fullname << endl;
    cout << " (" << abbrev2.Coords.Lat << ", " << abbrev2.Coords.Lon << ")" << endl;
    //Find the midpoint and closest building to the midpoint
    Coordinates middle = centerBetween2Points(abbrev1.Coords.Lat, abbrev1.Coords.Lon, abbrev2.Coords.Lat, abbrev2.Coords.Lon);
    BuildingInfo closestBuilding;
    double minimumn = numeric_limits<double>::max();
    for (auto& building : Buildings) {
      double currentDistance = distBetween2Points(middle.Lat, middle.Lon, building.Coords.Lat, building.Coords.Lon);
      if (currentDistance < minimumn) {
        minimumn = currentDistance;
        closestBuilding = building;
      }
    }
    BuildingInfo nearestBuilding;
    //Finds the closestBuilding to the current node
    nearestBuilding = findClosetBuilding(middle, Buildings);
    long long firstNode = -1;
    long long secondNode = -1;
    long long thirdNode = -1;
    //Finds the nearest nodes to each person and the destination building
    double minimum = INFINITY_VAL;
    double minimum1 = INFINITY_VAL;
    double minimum2 = INFINITY_VAL;

    cout << "Destination Building:" << endl;
    cout << " " << nearestBuilding.Fullname << endl;
    cout << " (" << nearestBuilding.Coords.Lat << ", " << nearestBuilding.Coords.Lon << ")" << endl << endl;
    //Finds the nearest node to person 1
    for (auto& vectFoot : Footways){ // Iterate over each footway in the vector
      for (int i = 0; i < vectFoot.Nodes.size(); i++){
        long long nodeOne = vectFoot.Nodes[i]; // Get the ID of the current node 
        Coordinates one = Nodes.at(nodeOne);  // Get the coordinates of the current node using the node ID
        double distanceOne = distBetween2Points(one.Lat, one.Lon, abbrev1.Coords.Lat, abbrev1.Coords.Lon);
        if (distanceOne < minimum) {  // Check if the calculated distance is less than the current minimum distance
          minimum = distanceOne; //Updates the minimum distance and the corresponding node ID
          firstNode = one.ID; 
        }
      }
    }
    cout << "Nearest P1 node:" << endl << " " << firstNode << endl;
    cout << " (" << Nodes[firstNode].Lat << ", " << Nodes[firstNode].Lon << ")" << endl;
    for (auto& vectFoot : Footways){ //Finds the nearest node to person 2
      for (int i = 0; i < vectFoot.Nodes.size(); i++){
        long long nodeTwo = vectFoot.Nodes[i];
        Coordinates two = Nodes.at(nodeTwo);
        // Calculate the distance between the current node and person 2's building
        double distanceTwo = distBetween2Points(two.Lat, two.Lon, abbrev2.Coords.Lat, abbrev2.Coords.Lon);
        if (distanceTwo < minimum1){
          minimum1 = distanceTwo;
          secondNode = two.ID;
        }
      }
    }
    cout << "Nearest P2 node:" << endl << " " << secondNode << endl;
    cout << " (" << Nodes[secondNode].Lat << ", " << Nodes[secondNode].Lon << ")" << endl;
    //finds the nearest node to the destination building
    for (auto& vectFoot : Footways) {
      for (int i = 0; i < vectFoot.Nodes.size(); i++){
        long long nodeThree = vectFoot.Nodes[i];
        Coordinates three = Nodes.at(nodeThree); // Get the coordinates of the current node using the node ID
        // Calculate the distance between the current node and the destination building
        double distanceThree = distBetween2Points(three.Lat, three.Lon, nearestBuilding.Coords.Lat, nearestBuilding.Coords.Lon);
        if (distanceThree < minimum2){
          minimum2 = distanceThree;
          thirdNode = three.ID;// update the minimum distance and the corresponding node ID
        }
      }
    }
    cout << "Nearest destination node:" << endl;
    cout << " " << thirdNode << endl;
    cout << " (" << Nodes[thirdNode].Lat << ", " << Nodes[thirdNode].Lon << ")" << endl;
    // Finds the shortest path and distances for person 1 and person 2 to the destination
    FindShortestPaths(G, firstNode, predOne, distanceOne);
    FindShortestPaths(G, secondNode, predTwo, distanceTwo);
    double distanceP1 = distanceOne[thirdNode]; // Gets the distances from the starting nodes to the destination node
    double distanceP2 = distanceTwo[thirdNode];

    if (distanceOne[secondNode] >= INFINITY_VAL) { // Check if there is no path from the starting node to the destination
      cout << "Sorry, destination unreachable." << endl;
    }
    if ((distanceP1 != INFINITY_VAL && distanceP2 != INFINITY_VAL)) { // Check if both Person 1 and Person 2 have a valid path to the destination
      cout << "Person 1's distance to dest " << distanceP1 << " miles" << endl;
      cout << "Path: ";
      stack <long long> dots; // Create a stack to store the nodes in the path
      dots.push(thirdNode); // Push the destination node onto the stack
      long long center = thirdNode; // Initialize the current node to the destination node
      while (predOne[center] != -1) { // Iterate through the predecessors to build the path on the stack
        dots.push(predOne[center]);
        center = predOne[center];
      }
      while (!dots.empty()){ // Output the path for Person 1
        center = dots.top();
        dots.pop();
        if (dots.size() >= 1){
          cout << center << "->";  // Output node followed by arrow if there are more nodes in the path
        }
        else {
          cout << center; // Output the last node in the path
          cout << endl;
        }
      }  
      cout << "Person 2's distance to dest " << distanceP2 << " miles" << endl;
      cout << "Path: ";
      stack <long long> dotsOne; // Create a stack to store the nodes in the path
      dotsOne.push(thirdNode);
      while (predTwo[center] != -1){
        dotsOne.push(predTwo[center]);
        center = predTwo[center];
      }
      while (!dotsOne.empty()){ // Output the path for Person 2
        center = dotsOne.top();
        dotsOne.pop(); //Removes the nodes off from the stack
        if(dotsOne.size() >= 1){
          cout << center << "->";
        }
        else {
          cout << center;
        }
      }
      cout << endl;
    }
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

int main() {
  graph<long long, double> G;

  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  for (const auto& node : Nodes) { // Add vertices to the graph for each node in the map Nodes
    G.addVertex(node.first);
  }

  for (const auto& footway : Footways) { // Add edges to the graph based on the coordinates of consecutive nodes in each footway
    const auto& nodes = footway.Nodes;
    for (size_t i = 0; i < nodes.size() - 1; ++i) {
      long long location1 = nodes[i];
      long long location2 = nodes[i + 1];
      // Calculate the distance between consecutive nodes
      double difference = distBetween2Points(Nodes[location1].Lat, Nodes[location1].Lon, Nodes[location2].Lat, Nodes[location2].Lon);
      // Add edges to the graph with the calculated distance
      G.addEdge(location1, location2, difference);
      G.addEdge(location2, location1, difference);
    }
  }

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  // Execute Application
  application(Nodes, Footways, Buildings, G);

  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}