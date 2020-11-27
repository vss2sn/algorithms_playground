#ifndef GRAPHAM_HPP
#define GRAPHAM_HPP

#include <algorithm>
#include <random>
#include <stack>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "algorithms/graph.hpp"

namespace graphAM {

// Double is to allow non integral weights
// using GraphAM = std::vector<std::vector<double>>;
// using GraphAL = std::vector<std::vector<std::pair<int, double>>>;
// using GraphEL = std::vector<Edge>;

class GraphAM : public Graph {
private:
  std::vector<std::vector<double>> g;
  int v;

  template <typename T>
  void GraphAMUtil(const int V, T& dist);

  template <typename T>
  void fill_row(std::vector<double>& row, T& dist);

  double fordFulkersonRGUtil(int source, int sink);

  int UnionFindFindUtil(const int v, const std::vector<int>& parent) const;
  void UnionFindUnionUtil(const int v1, const int v2, std::vector<int>& parent) const;

  int UnionFindRCFindUtil(const int v, std::vector<std::pair<int, int>>& subsets) const;
  void UnionFindRCUnionUtil(const int v1, const int v2, std::vector<std::pair<int, int>>& subsets) const;

  void DFSUtil(const int source, std::vector<bool>& visited) const;

  bool colourGraphUtil(const int vert, const int n_c, std::vector<int>& colour) const;
  void allPathsBetweenUtil(const int source, int sink, std::vector<int>& path,
      std::vector<std::vector<int>>& paths, std::vector<bool>& visited) const;

  void PrintAllDijkstraPathsFound(const int source, const std::vector<std::pair<int, double>>& p_md) const;
  void DFSUtilWithFinishTime(int source, std::vector<bool>& visited, std::stack<int>& visit_order) const;

  std::vector<std::vector<int>> KosarajuAlgorithmUtil();
  bool StronglyConnectedKosarajuUtil();

  std::vector<int> HierholzersAlgorithmUtil();

  void ArticulationPointsUtil(const int vert1, std::vector<bool>& visited, std::vector<int>& parent,
    std::vector<int>& tod, std::vector<int>& low, std::unordered_set<int>& articulation_points, int time) const ;

  bool TopologicalSortUtil(std::vector<int>& sorted,
    std::vector<bool>& visited, std::vector<bool>& this_cycle,
    const int vert) const;

  void FindBridgesUtil(const int vert1, std::vector<bool>& visited, std::vector<int>& parent,
    std::vector<int>& tod, std::vector<int>& low, std::vector<std::pair<int, int>>& bridges, int time) const;
public:
  GraphAM(const int V = 11, const bool random = true, const bool different_weights = true);
  GraphAM(const GraphAM& g_am);
  template <typename T>
  GraphAM(const int V, T& dist);
  GraphAM(const std::vector<std::vector<double>>& g_v);

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

  // Find mother Vertex
  virtual std::tuple<bool, int> FindMotherVertex () const override;

  // All Paths Between 2 vertices
  virtual std::tuple<bool, std::vector<std::vector<int>>> allPathsBetween(const int source, const int sink) const override;

  // Colour graph using n_c colours
  virtual std::tuple<bool, std::vector<int>> colourGraph(int n_c) const override;

  // Dijkstra from source to sink
  virtual std::tuple<bool, double, std::vector<int>> Dijkstra(const int source, const int sink) const override;

  // Dijkstra from source to all possible points
  virtual std::tuple<bool, std::vector<std::pair<int, double>>> Dijkstra(const int source) const override;

  // Find all strongly connected components of a graph
  virtual std::vector<std::vector<int>> KosarajuAlgorithm() const override;

  // Check if graph is strongly connected
  virtual bool StronglyConnectedKosaraju() const override;

  // Creates minimum spanning tree for graph
  virtual std::vector<Edge> KruskalsAlgorithm() const override;

  // Find minimum distance between every pair of points
  virtual std::vector<std::vector<double>> FloydWarshall() const override;

  // Detect negative cycle
  virtual bool BellmanFord(const int source) const override;

  // Find Eulerian Path
  virtual std::vector<int> HierholzersAlgorithm() const override;

  // Find articulation points in an undirected graph
  virtual std::unordered_set<int> ArticulationPoints() const override;

  // Topological sort (checks whether the graph is a directed acyclic as well)
  virtual std::tuple<bool, std::vector<int>> TopologicalSort() const override;

  virtual std::vector<std::pair<int, int>> FindBridges() const override;

  virtual bool IsBipartite() const override;

  virtual bool DivideIntoTwoCliques() const override;

};

/* ------------------------------- */
/* Templated functions definitions */
/* ------------------------------- */

template<typename T>
void GraphAM::fill_row(std::vector<double>& row, T& dist)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::generate(row.begin(), row.end(), [&]() {return dist(gen);});
}

template <typename T>
void GraphAM::GraphAMUtil(const int V, T& dist) {
  v = V;
  g = std::vector<std::vector<double>>(V, std::vector<double>(V, 0));
  std::for_each(g.begin(), g.end(), [&](std::vector<double>& row){ fill_row(row, dist); });
}

template <typename T>
GraphAM::GraphAM(const int V, T& dist) {
  GraphAMUtil(V, dist);
}

}  // namespace graphAM

#endif  // GRAPHAM_HPP
