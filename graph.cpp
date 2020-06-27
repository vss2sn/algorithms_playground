#include <iostream>

#include "graph.hpp"

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
bool Graph::Prim() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return false;
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
