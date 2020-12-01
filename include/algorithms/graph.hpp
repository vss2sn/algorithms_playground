#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <tuple>
#include <unordered_set>
#include <vector>

/**
* @brief Class Edge
* @details Contains the start vertex, end vertex, and the weight of the edge
*/
struct Edge {
  int u, v;
  double w;

  bool operator == (const Edge& e2) const;
};

/**
* @brief Class Graph
* @details Abstract class from which other graph class are derived
*/
class Graph {
public:

  /**
   * @brief Prints the graph
   */
  virtual void PrintGraph() const;

  /**
   * @brief Breadth first search
   * @param [in] source source vertex
   * @param [in] sink sink vertex
   * @return bool whether the path exists and the path itself
   */
  virtual std::tuple<bool, std::vector<int>> BFS (int source , int sink) const;

  /**
   * @brief Depth first search
   * @param [in] source source vertex
   * @param [in] sink sink vertex
   * @return bool whether the path exists and the path itself
   */
  virtual std::tuple<bool, std::vector<int>> DFS (int source , int sink) const;

  /**
   * @brief Prim's algorithm for finding the minimum spanning tree
   * @return bool for whether the algorithm was successful and the minimum spanning tree
   * @details Prim's algorithm only works on undirected complete graphs
   */
  virtual std::tuple<bool, std::vector<Edge>> Prim() const;

  /**
   * @brief Ford Fulkerson algorithm for finding the maximum flow through a graph
   * @param [in] source vertex from which the flows originate
   * @param [in] sink vertex into which the flows flow
   * @return maximum flow value
   */
  virtual double fordFulkerson(int source, int sink) const;

  /**
   * @brief Detect cycles within a graph using the UnionFind algorithm
   * @return bool whether a cycle is detected
   */
  virtual bool UnionFindDetectCycle() const;

  /**
   * @brief Detect cycles within a graph using the UnionFind algorithm with Rank Compression
   * @return bool whether a cycle is detected
   */
  virtual bool UnionFindRCDetectCycle() const;

  /**
   * @brief Find the mother vertex of a graph
   * @return bool whether a mother vertex was detected and said vertex
   */
  virtual std::tuple<bool, int> FindMotherVertex() const;

  /**
   * @brief Finds all paths between a source and a sink
   * @param [in] source source vertex
   * @param [in] sink sink vertex
   * @return bool whether paths exit and a vector of all the paths between the source and the sink
   */
  virtual std::tuple<bool, std::vector<std::vector<int>>> allPathsBetween(const int source, const int sink) const;

  /**
   * @brief Colour a graph
   * @param [in] n_c number of colours that can be used to colour the graph
   * @return bool whether the gaph was coloured and a vector representing the colour of each vertex
   */
  virtual std::tuple<bool, std::vector<int>> colourGraph(int n_c) const;

  /**
   * @brief Dijkstra's algorithm
   * @param [in] source source vertex
   * @param [in] sink sink vertex
   * @return bool whether a path was found and said path
   */
  virtual std::tuple<bool, double, std::vector<int>> Dijkstra(const int source, const int sink) const;

  /**
   * @brief Dijkstra's algorithm
   * @param [in] source source vertex
   * @return bool whether paths were found and said paths
   */
  virtual std::tuple<bool, std::vector<std::pair<int, double>>> Dijkstra(const int source) const;

  /**
   * @brief Kosaraju's algorithm
   * @return vector of strongly connnected components
   */
  virtual std::vector<std::vector<int>> KosarajuAlgorithm() const;

  /**
   * @brief Check if graph is strongly connected using Kosaraju's algorithm
   * @return bool whether the graph is strongly connected
   */
  virtual bool StronglyConnectedKosaraju() const;

  /**
   * @brief Kruskal's algorithm for finding the minimum spanning tree
   * @return vector of edges within the minimum spanning tree
   */
  virtual std::vector<Edge> KruskalsAlgorithm() const;

  /**
   * @brief Floyd Warshall's algorithm to find the minimum distance between every pair of points
   * @return vector of paths
   */
  virtual std::vector<std::vector<double>> FloydWarshall() const;

  /**
   * @brief Detect negative cycle within the graph
   * @param [in] source source vertex
   * @return bool whether the graph has a negative cycle
   */
  virtual bool BellmanFord(const int source) const;

  /**
   * @brief Hierholzer's algorithm find eulerian paths
   * @return returns an eulerian path
   */
  virtual std::vector<int> HierholzersAlgorithm() const;

  /**
   * @brief Find articulation points in an undirected graph
   * @return unordered set of articulation points
   */
  virtual std::unordered_set<int> ArticulationPoints() const;

  /**
   * @brief Find a topological sorted order of vertices of the graph
   * @return bool whether a topological sorted order exists and a topological sorted ordering of the vertices
   * @details checks whether the graph is a directed acyclic as well
   */
  virtual std::tuple<bool, std::vector<int>> TopologicalSort() const;

  /**
   * @brief Find bridges within the graph
   * @return vector of edges represented as start and end vertices that are bridges
   */
  virtual std::vector<std::pair<int, int>> FindBridges() const;

  /**
   * @brief checks whether the graph is bipartite
   * @return bool whether the  graph is bipartite
   */
  virtual bool IsBipartite() const;

  /**
   * @brief checks whether the graph can be divided into 2 cliques
   * @return bool whether the graph can be divided into 2 cliques
   */
  virtual bool DivideIntoTwoCliques() const;

  /**
   * @brief Assigns each vertex a level based on how many edges it took to reach the vertex from the source and checks whether the sink was reached
   * @param [in] source source vertex
   * @param [in] sink sink vertex
   * @return bool whether the graph levels were assigned and the level of each vertex
   */
  virtual std::tuple<bool, std::vector<int>> CreateLevelGraph(const int source, const int sink) const;

  /**
   * @brief Assigns each vertex a level based on how many edges it took to reach the vertex from the source and checks whether the sink was reached
   * @param [in] source source vertex
   * @param [in] sink sink vertex
   * @param [in] g graph for which the levels are to be assigned
   * @return bool whether the graph levels were assigned and the level of each vertex
   */
  virtual std::tuple<bool, std::vector<int>> CreateLevelGraph(const int source, const int sink, const std::vector<std::vector<double>>& g) const;

  /**
   * @brief Calculate maximum flow between a source and a sink using Dinac's algorithm
   * @param [in] source source vertex
   * @param [in] sink sink vertex
   * @return maximum flow
   */
  virtual double DinacsAlgorithm(const int source, const int sink) const;

  /**
   * @brief Find a Hamiltonian path if one exists
   * @return bool whether a Hamniltonian path was found and the path
   */
  virtual std::tuple<bool, std::vector<int>> HamiltonianPath() const;
};

#endif  // GRAPH_HPP
