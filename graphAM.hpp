#ifndef GRAPHAM_HPP
#define GRAPHAM_HPP

#include <algorithm>
#include <random>
#include <tuple>
#include <vector>

#include "graph.hpp"

namespace graphAM {

struct Edge {
  int u, v;
  double w;

  Edge(int u, int v, double w);
  bool operator < (const Edge& other) const;
  bool operator <= (const Edge& other) const;
  bool operator > (const Edge& other) const;
  bool operator >= (const Edge& other) const;

};

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
  void UnionFindUnionUtil(int v1, int v2, std::vector<int>& parent) const;

  int UnionFindRCFindUtil(int v, std::vector<std::pair<int, int>>& subsets) const;
  void UnionFindRCUnionUtil(int v1, int v2, std::vector<std::pair<int, int>>& subsets) const;

  void DFSUtil(const int source, std::vector<bool>& visited) const;

  bool colourGraphUtil(const int vert, const int n_c, std::vector<int>& colour) const;
  void allPathsBetweenUtil(const int source, int sink, std::vector<int>& path,
      std::vector<std::vector<int>>& paths, std::vector<bool>& visited) const;

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
  virtual bool Prim() const override;

  // Ford Fulkerson
  virtual double fordFulkerson(int source, int sink) const override;

  // Union Find Algorithms
  virtual bool UnionFindDetectCycle() const override;
  virtual bool UnionFindRCDetectCycle() const override;

  // Find mother Vertex
  virtual std::tuple<bool, int> FindMotherVertex () const override;

  // All Paths Between 2 vertices
  virtual std::tuple<bool, std::vector<std::vector<int>>> allPathsBetween(int source, int sink) const override;

  // Colour graph using n_c colours
  virtual std::tuple<bool, std::vector<int>> colourGraph(int n_c) const override;


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