#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <tuple>
#include <unordered_set>
#include <vector>

struct Edge {
  int u, v;
  double w;

  bool operator == (const Edge& e2) const;
};

class Graph {
public:

  virtual void PrintGraph() const;
  virtual std::tuple<bool, std::vector<int>> BFS (int source , int sink) const;
  virtual std::tuple<bool, std::vector<int>> DFS (int source , int sink) const;

  // Prim's algorithm only works on undirected complete graphs
  virtual std::tuple<bool, std::vector<Edge>> Prim() const;

  // Ford Fulkerson
  virtual double fordFulkerson(int source, int sink) const;

  // Union Find Algorithms
  virtual bool UnionFindDetectCycle() const;
  virtual bool UnionFindRCDetectCycle() const;

  // Find mother Vertex
  virtual std::tuple<bool, int> FindMotherVertex () const;

  // All Paths Between 2 vertices
  virtual std::tuple<bool, std::vector<std::vector<int>>> allPathsBetween(int source, int sink) const;

  // Colour graph using n_c colours
  virtual std::tuple<bool, std::vector<int>> colourGraph(int n_c) const;

  // Dijkstra from source to sink
  virtual std::tuple<bool, double, std::vector<int>> Dijkstra(const int source, const int sink) const;

  // Dijkstra from source to all possible points
  virtual std::tuple<bool, std::vector<std::pair<int, double>>> Dijkstra(const int source) const;

  // Find all strongly connected components of a graph
  virtual std::vector<std::vector<int>> KosarajuAlgorithm() const;

  // Check if graph is strongly connected
  virtual bool StronglyConnectedKosaraju() const;

  // Creates minimum spanning tree for graph
  virtual std::vector<Edge> KruskalsAlgorithm() const;

  // Find minimum distance between every pair of points
  virtual std::vector<std::vector<double>> FloydWarshall() const;

  // Detect negative cycle
  virtual bool BellmanFord(const int source) const;

  // Find Eulerian Paths
  virtual std::vector<int> HierholzersAlgorithm() const;

  // Find articulation points in an undirected graph
  virtual std::unordered_set<int> ArticulationPoints() const;

  // Topological sort (checks whether the graph is a directed acyclic as well)
  virtual std::tuple<bool, std::vector<int>> TopologicalSort() const;
};

#endif  // GRAPH_HPP
