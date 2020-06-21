#ifndef GRAPHAM_HPP
#define GRAPHAM_HPP

#include <algorithm>
#include <random>
#include <tuple>
#include <vector>

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
using GraphAM = std::vector<std::vector<double>>;
using GraphAL = std::vector<std::vector<std::pair<int, double>>>;
using GraphEL = std::vector<Edge>;

void PrintGraph(const GraphAM& g);
GraphAM createGraphAM(const int V = 11, const bool random = true, const bool different_weights = true);
std::tuple<bool, std::vector<int>> BFS (GraphAM& g, int source , int sink);
std::tuple<bool, std::vector<int>> DFS (GraphAM& g, int source , int sink);

// Prim's algorithm only works on undirected complete graphs
bool Prim(GraphAM& g);

// Ford Fulkerson
double fordFulkerson(GraphAM& g, int source, int sink);

// Union Find Algorithms
bool UnionFindDetectCycle(const GraphAM& g);
bool UnionFindRCDetectCycle(const GraphAM& g);

// Find mother Vertex
std::tuple<bool, int> FindMotherVertex (const GraphAM& g);

// All Paths Between 2 vertices
std::tuple<bool, std::vector<std::vector<int>>> allPathsBetween(const GraphAM& g, int source, int sink);

// Colour graph using n_c colours
std::tuple<bool, std::vector<int>> colourGraph(const GraphAM& g, int n_c);

/* ------------------------------- */
/* Templated functions definitions */
/* ------------------------------- */

template<typename T>
void fill_row(std::vector<double>& row, T& dist)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::generate(row.begin(), row.end(), [&]() {return dist(gen);});
}

template <typename T>
GraphAM createGraphAMUtil(const int V, T& dist) {
  GraphAM g(V, std::vector<double>(V, 0));
  std::for_each(g.begin(), g.end(), [&](std::vector<double>& row){ fill_row(row, dist); });
  return g;
}

template <typename T>
GraphAM createGraphAM(const int V, T& dist) {
  return createGraphAMUtil(V, dist);
}

}  // namespace graphAM

#endif  // GRAPHAM_HPP
