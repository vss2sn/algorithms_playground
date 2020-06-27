#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <tuple>
#include <vector>

class Graph {
public:

  virtual void PrintGraph() const;
  virtual std::tuple<bool, std::vector<int>> BFS (int source , int sink) const;
  virtual std::tuple<bool, std::vector<int>> DFS (int source , int sink) const;

  // Prim's algorithm only works on undirected complete graphs
  virtual bool Prim() const;

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

};

#endif  // GRAPH_HPP
