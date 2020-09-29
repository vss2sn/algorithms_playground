#include <iostream>

#include "algorithms/graph.hpp"

bool Edge::operator == (const Edge& e2) const {
  return (this->u == e2.u && this->v == e2.v && this->w == e2.w);
}

void Graph::PrintGraph() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
}

std::tuple<bool, std::vector<int>> Graph::BFS (int source , int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, std::vector<int>()};
}
std::tuple<bool, std::vector<int>> Graph::DFS (int source , int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, std::vector<int>()};
}

// Prim's algorithm only works on undirected complete graphs
std::tuple<bool, std::vector<Edge>> Graph::Prim() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, std::vector<Edge>()};
}

// Ford Fulkerson
double Graph::fordFulkerson(int source, int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return 0;
}

// Union Find Algorithms
bool Graph::UnionFindDetectCycle() const{
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return false;
}

bool Graph::UnionFindRCDetectCycle() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return false;
}

// Find mother Vertex
std::tuple<bool, int> Graph::FindMotherVertex () const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, 0};
}

// All Paths Between 2 vertices
std::tuple<bool, std::vector<std::vector<int>>> Graph::allPathsBetween(int source, int sink) const{
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, std::vector<std::vector<int>>()};
}

// Colour graph using n_c colours
std::tuple<bool, std::vector<int>> Graph::colourGraph(int n_c) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, std::vector<int>()};
}

std::tuple<bool, double, std::vector<int>> Graph::Dijkstra(const int source, const int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, 0, std::vector<int>()};
}

// Dijkstra from source to all possible points
std::tuple<bool, std::vector<std::pair<int, double>>> Graph::Dijkstra(const int source) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, std::vector<std::pair<int, double>>()};
}

// Find all strongly connected components of a graph
std::vector<std::vector<int>> Graph::KosarajuAlgorithm() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<std::vector<int>>();
}

// Check if graph is strongly connected
bool Graph::StronglyConnectedKosaraju() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return false;
}

// Creates minimum spanning tree for graph
std::vector<Edge> Graph::KruskalsAlgorithm() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<Edge>();
}

// Find minimum distance between every pair of points
std::vector<std::vector<double>> Graph::FloydWarshall() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<std::vector<double>>();
 }

 // Detect negative cycle
 bool Graph::BellmanFord(const int source) const {
   std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return false;
}
