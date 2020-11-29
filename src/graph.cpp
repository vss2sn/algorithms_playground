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

std::tuple<bool, std::vector<Edge>> Graph::Prim() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, std::vector<Edge>()};
}

double Graph::fordFulkerson(int source, int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return 0;
}

bool Graph::UnionFindDetectCycle() const{
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return false;
}

bool Graph::UnionFindRCDetectCycle() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return false;
}

std::tuple<bool, int> Graph::FindMotherVertex () const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, 0};
}

std::tuple<bool, std::vector<std::vector<int>>> Graph::allPathsBetween(const int source, const int sink) const{
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, std::vector<std::vector<int>>()};
}

std::tuple<bool, std::vector<int>> Graph::colourGraph(int n_c) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, std::vector<int>()};
}

std::tuple<bool, double, std::vector<int>> Graph::Dijkstra(const int source, const int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, 0, std::vector<int>()};
}

std::tuple<bool, std::vector<std::pair<int, double>>> Graph::Dijkstra(const int source) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return {0, std::vector<std::pair<int, double>>()};
}

std::vector<std::vector<int>> Graph::KosarajuAlgorithm() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<std::vector<int>>();
}

bool Graph::StronglyConnectedKosaraju() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return false;
}

std::vector<Edge> Graph::KruskalsAlgorithm() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<Edge>();
}

std::vector<std::vector<double>> Graph::FloydWarshall() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<std::vector<double>>();
 }

 bool Graph::BellmanFord(const int source) const {
   std::cout << __FUNCTION__ << " not yet defined" << '\n';
   return false;
}

std::vector<int> Graph::HierholzersAlgorithm() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<int>();
}

std::unordered_set<int> Graph::ArticulationPoints() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::unordered_set<int>();
}

std::tuple<bool, std::vector<int>> Graph::TopologicalSort() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, std::vector<int>()};
}

std::vector<std::pair<int, int>> Graph::FindBridges() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return std::vector<std::pair<int, int>>();
}

bool Graph::IsBipartite() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return false;
}

bool Graph::DivideIntoTwoCliques() const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return false;
}

std::tuple<bool, std::vector<int>> Graph::CreateLevelGraph(const int source, const int sink, const std::vector<std::vector<double>>& g) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, std::vector<int>()};
}

std::tuple<bool, std::vector<int>> Graph::CreateLevelGraph(const int source, const int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return {false, std::vector<int>()};
}

double Graph::DinacsAlgorithm(const int source, const int sink) const {
  std::cout << __FUNCTION__ << " not yet defined" << '\n';
  return 0;
}
