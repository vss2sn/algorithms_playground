#ifndef GRAPHAL_HPP
#define GRAPHAL_HPP

#include <algorithm>
#include <random>
#include <stack>
#include <tuple>
#include <vector>
#include <utility>

#include "graph.hpp"

namespace graphAL {

class GraphAL : public Graph {
private:
  std::vector<std::vector<std::pair<int, double>>> g;
  int v;

  // template <typename T>
  // void GraphALUtil(const int V, T& dist);
  //
  // template <typename T>
  // void fill_row(std::vector<double>& row, T& dist);
  
  double fordFulkersonRGUtil(int source, int sink);

  int UnionFindFindUtil(const int v, const std::vector<int>& parent) const;
  void UnionFindUnionUtil(int v1, int v2, std::vector<int>& parent) const;

  int UnionFindRCFindUtil(int v, std::vector<std::pair<int, int>>& subsets) const;
  void UnionFindRCUnionUtil(int v1, int v2, std::vector<std::pair<int, int>>& subsets) const;

  // void DFSUtil(const int source, std::vector<bool>& visited) const;
  //
  // bool colourGraphUtil(const int vert, const int n_c, std::vector<int>& colour) const;
  // void allPathsBetweenUtil(const int source, int sink, std::vector<int>& path,
  //     std::vector<std::vector<int>>& paths, std::vector<bool>& visited) const;
  //
  // void PrintAllDijkstraPathsFound(const int source, const std::vector<std::pair<int, double>>& p_md) const;
  // void DFSUtilWithFinishTime(int source, std::vector<bool>& visited, std::stack<int>& visit_order) const;
  //
  // std::vector<std::vector<int>> KosarajuAlgorithmUtil();
  // bool StronglyConnectedKosarajuUtil();

public:
  // GraphAL(const int V = 11, const bool random = true, const bool different_weights = true);
  // GraphAL(const GraphAL& g_am);
  // template <typename T>
  // GraphAL(const int V, T& dist);
  GraphAL(const std::vector<std::vector<std::pair<int, double>>>& g_v);

  virtual void PrintGraph() const override;
  virtual std::tuple<bool, std::vector<int>> BFS (int source , int sink) const override;
  virtual std::tuple<bool, std::vector<int>> DFS (int source , int sink) const override;

  // Prim's algorithm only works on undirected complete graphs
  virtual std::tuple<bool, std::vector<Edge>> Prim() const override;

  // Ford Fulkerson
  virtual double fordFulkerson(int source, int sink) const override;

  // Union Find Algorithms
  virtual bool UnionFindDetectCycle() const override;
  virtual bool UnionFindRCDetectCycle() const override;
//
//   // Find mother Vertex
//   virtual std::tuple<bool, int> FindMotherVertex () const override;
//
//   // All Paths Between 2 vertices
//   virtual std::tuple<bool, std::vector<std::vector<int>>> allPathsBetween(int source, int sink) const override;
//
//   // Colour graph using n_c colours
//   virtual std::tuple<bool, std::vector<int>> colourGraph(int n_c) const override;
//
//   // Dijkstra from source to sink
//   virtual std::tuple<bool, double, std::vector<int>> Dijkstra(const int source, const int sink) const override;
//
//   // Dijkstra from source to all possible points
//   virtual std::tuple<bool, std::vector<std::pair<int, double>>> Dijkstra(const int source) const override;
//
//   // Find all strongly connected components of a graph
//   virtual std::vector<std::vector<int>> KosarajuAlgorithm() const override;
//
//   // Check if graph is strongly connected
//   virtual bool StronglyConnectedKosaraju() const override;
//
  // Creates minimum spanning tree for graph
  virtual std::vector<Edge> KruskalsAlgorithm() const override;
//
//   // Find minimum distance between every pair of points
//   virtual std::vector<std::vector<double>> FloydWarshall() const override;
//
//   // Detect negative cycle
//   virtual bool BellmanFord(const int source) const override;
};

}  // namespace graphAL

#endif  // GRAPHAL_HPP
