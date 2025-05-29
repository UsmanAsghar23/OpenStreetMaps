// graph.h
#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <iterator>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
     map<VertexT, map<VertexT, WeightT>> adjacencyList;

  //
  // We are using adjacency matrix implementation, where rows
  // are the starting vertex and cols are the ending vertex.
  // We keep track of the vertices in the Vertices vector,
  // where the vertex's position in the vector --- 0, 1, 2,
  // 3, 4, 5, ... --- denotes the row in the adjacency matrix
  // where their edges are found.  Example: if vertex "ORD" is
  // in position 1 of the Vertices vector, then row 1 of
  // AdjMatrix are the edges that start at "ORD" and lead to
  // other vertices.
  //
  // static constexpr int MatrixSize = 100;

  // EdgeData         AdjMatrix[MatrixSize][MatrixSize];
  // vector<VertexT>  Vertices;

  //
  // _LookupVertex
  //
  // Finds the vertex in the Vertices vector and returns it's
  // index position if found, otherwise returns -1.
  //
  // int _LookupVertex(VertexT v) const {
  //   for (int i = 0; i < this->NumVertices(); ++i) {
  //     if (this->Vertices[i] == v)  // already in the graph:
  //       return i;
  //   }

    // if get here, not found:
    // return -1;
  // }

 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //
  graph() {

  }

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const {
    return adjacencyList.size(); 
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    int count = 0;
    for (const auto& pair : this->adjacencyList) { 
      count += (pair.second.size());
    }
    return count;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    // Check if the vertex already exists
    if (this->adjacencyList.find(v) != this->adjacencyList.end()) {
      return false;
    }
    // Add the vertex to the graph
    this->adjacencyList[v] = map<VertexT, WeightT>();
    return true;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    // Check if the vertices exist
    if (this->adjacencyList.find(from) == this->adjacencyList.end() || this->adjacencyList.find(to) == this->adjacencyList.end()) {
      return false;
    }
    // Add the edge to the adjacency list
    this->adjacencyList[from][to] = weight;
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    // Check if the vertices and edge exist
    if (this->adjacencyList.find(from) == this->adjacencyList.end() || this->adjacencyList.find(to) == this->adjacencyList.end()) {
      return false;
    }
    // Find the edge and return the weight
    auto it = this->adjacencyList.at(from).find(to);
    if (it != this->adjacencyList.at(from).end()) {
      weight = it->second;
      return true;
    }
    return false;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;
    // Check if the vertex exists
    if (this->adjacencyList.find(v) == this->adjacencyList.end()) {
      return S;
    }
    // Add neighbors to the set
    for (const auto& pair : this->adjacencyList.at(v)) {
      S.insert(pair.first);
    }
    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector<VertexT> vertices;
    for (const auto& pair : adjacencyList) {
      vertices.push_back(pair.first);
    }
    return vertices;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;
    output << "**Num vertices: " << NumVertices() << endl;
    output << "**Num edges: " << NumEdges() << endl;
    output << endl;
    output << "**Vertices:" << endl;

    for (const auto& vertex : adjacencyList) {
      output << " " <<vertex.first << endl;
    }
    output << endl;
    output << "**Edges:" << endl;

    for (const auto& pair : adjacencyList) {
      output << " " << pair.first << ": ";
      for (const auto& edge : pair.second) {
        output << "(" << edge.first << ", " << edge.second << ") ";
      }
      output << endl;
    }
    output << "**************************************************" << endl;
  }

  void clear() {
    this->adjacencyList.clear();
  }
};
