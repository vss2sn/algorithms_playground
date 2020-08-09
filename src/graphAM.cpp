#include <algorithm>
#include <iostream>
#include <queue>
#include <limits>
#include <random>
#include <stack>
#include <tuple>
#include <utility>
#include <vector>

#include "algorithms/graphAM.hpp"

namespace graphAM {

void GraphAM::PrintGraph() const {
  std::cout << "GraphAM: " << '\n';
  for(const auto& row : g) {
    for(const auto ele : row) {
      std::cout << ele << ' ';
    }
    std::cout << '\n';
  }
}

GraphAM::GraphAM(const GraphAM& g_am) {
  this->g = g_am.g;
  this->v = g_am.v;
}
GraphAM::GraphAM(const int V, const bool random, const bool different_weights) {
  g = std::vector<std::vector<double>>(V, std::vector<double>(V,0));
  if(random) {
    if(different_weights){
      std::uniform_int_distribution<> dist(0, V);
      GraphAMUtil(V, dist);
    }
    else {
      std::uniform_int_distribution<> dist(0, 1);
      GraphAMUtil(V, dist);
    }
  }
  else {
    std::vector<std::vector<double>>(V, std::vector<double>(V, 0));
  }
}

GraphAM::GraphAM(const std::vector<std::vector<double>>& g_v) {
  g = g_v;
  v = g.size();
}


std::tuple<bool, std::vector<int>> GraphAM::BFS (int source , int sink) const {

  if(source >= v || sink >= v ) return {false, std::vector<int>()};

  std::vector<int> parent(v, -1);

  std::queue<int> q;
  q.push(source);
  parent[source] = source;

  while(!q.empty()) {
    if(parent[sink]!=-1) break;

    int current_vertex = q.front();
    q.pop();

    for(int r_index = 0; r_index < v; r_index++) {
      if(g[current_vertex][r_index] != 0 && parent[r_index] == -1) {
        parent[r_index] = current_vertex;
        if(r_index == sink) break;
        q.push(r_index);
      }
    }
  }

  if(parent[sink]!=-1) {
    std::vector<int> path;
    int current_vertex = sink;
    path.push_back(current_vertex);
    while (current_vertex!=source) {
      current_vertex = parent[current_vertex];
      path.push_back(current_vertex);
    }
    std::reverse(path.begin(), path.end());
    return {true, path};
  } else {
    return {false, std::vector<int>()};
  }

}

std::tuple<bool, std::vector<int>> GraphAM::DFS (int source , int sink) const {

  if(source >= v || sink >= v ) return {false, std::vector<int>()};

  std::vector<int> parent(v, -1);

  std::stack<int> s;
  s.push(source);
  parent[source] = source;

  while(!s.empty()) {
    if(parent[sink]!=-1) break;

    int current_vertex = s.top();
    s.pop();

    for(int r_index = 0; r_index < v; r_index++) {
      if(g[current_vertex][r_index] != 0 && parent[r_index] == -1) {
        parent[r_index] = current_vertex;
        if(r_index == sink) break;
        s.push(r_index);
      }
    }
  }

  if(parent[sink]!=-1) {
    std::vector<int> path;
    int current_vertex = sink;
    path.push_back(current_vertex);
    while (current_vertex!=source) {
      current_vertex = parent[current_vertex];
      path.push_back(current_vertex);
    }
    std::reverse(path.begin(), path.end());
    return {true, path};
  } else {
    return {false, std::vector<int>()};
  }

}

// Prim's algorithm only works on undirected complete graphs
// This will fail if the adjacency matrix is not symetrical
std::tuple<bool, std::vector<Edge>> GraphAM::Prim() const {
  if( v == 0 ) return {false, std::vector<Edge>()};

  for(int i =0; i < v; i++) {
    for (int j = i+1; j < v; j++) {
      if(g[i][j] != g[j][i]) {
        std::cout << "Edge from " << i << " to " << j << " not symetric/undirected. \
        Prim's algorithm  will fail. " << '\n';
        return {false, std::vector<Edge>()};
      }
    }
  }

  // This can be converted to a map of <point, parent>
  // BUt given that by the end of this algorithm it will contain v elements
  // anyway, sticking to vector for now
  // Vector size is unchanging so no copy operations; and access is O(1)
  std::vector<int> parent(v, -1);

  int source = 0;
  parent[source] = source;

  struct EdgeCompare {
    bool operator () (const Edge e1, const Edge e2) const {
      return e1.w > e2.w;
    }
  };

  std::priority_queue<Edge, std::vector<Edge>, EdgeCompare> edge_pq;

  // Insert current vertex
  for (int i = 0; i< v; i++) {
    if( g[source][i] != 0) {
      if (i!=source) edge_pq.push(Edge(source, i, g[source][i]));
    }
  }

  std::vector<Edge> mst;
  // Iterate until no edges left
  while(!edge_pq.empty()) {
    Edge e = edge_pq.top();
    edge_pq.pop();
    if(parent[e.v] == -1){
      parent[e.v] = e.u;
      mst.emplace_back(e);
      for (int i = 0; i< v; i++) {
        if( g[e.v][i] != 0 && parent[i] == -1 ) {
          edge_pq.push(Edge(e.v, i, g[e.v][i]));
        }
      }
    }
  }

  for(const auto& pid: parent) {
    if(pid ==-1 ) {
      std::cout << "Point " << pid << " not  included in MST " << '\n';
      return {false, std::vector<Edge>()};
    }
  }

  return {true, mst};
}

double GraphAM::fordFulkersonRGUtil(int source, int sink) {
  double max_flow = 0;
  while(true){
    auto [found_path, path] = BFS(source , sink);
    if(!found_path) break;
    double flow = std::numeric_limits<double>::max();
    for(int i = 1; i < path.size(); i++ ) {
      flow = std::min(flow, g[path[i-1]][path[i]]);
    }
    for(int i = 1; i < path.size(); i++ ) {
      g[path[i-1]][path[i]] -= flow;
      g[path[i]][path[i-1]] += flow;
    }
    max_flow +=flow;
  }
  return max_flow;
}

double GraphAM::fordFulkerson(int source, int sink) const {
  GraphAM rg(*this);
  return rg.fordFulkersonRGUtil(source, sink);
}


int GraphAM::UnionFindFindUtil(const int v, const std::vector<int>& parent) const {
  if(parent[v] == -1) return v;
  return UnionFindFindUtil(parent[v], parent);
}

void GraphAM::UnionFindUnionUtil(const int v1, const int v2, std::vector<int>& parent) const {
  // The two calls below to find are not required as in the union find algorithm,
  // the parents are already found and the unionutil is only called when they
  // are not equal. Left here for completeness (and future use) and as it's only going to be a
  // single call, not a recursive one as their parents will be -1
  int rv1 = UnionFindFindUtil(v1, parent);
  int rv2 = UnionFindFindUtil(v2, parent);
  if (rv1 != rv2) parent[rv2] = rv1;
}

bool GraphAM::UnionFindDetectCycle() const {

  std::vector<int> parent(v, -1);
  for (int i = 0; i < v ; i++) {
    int v1 = UnionFindFindUtil(i, parent);
    for (int j = 0; j < v; j++) {
      if (g[i][j] != 0) {
        int v2 = UnionFindFindUtil(j, parent);
        if(v1 == v2) {
          std::cout << __FUNCTION__ << " | " <<  " Cycle detected" << '\n';
          std::cout << __FUNCTION__ << " | " <<  " Edge vertices " << i << ' ' << j << '\n';
          std::cout << __FUNCTION__ << " | " <<  " Common vertex " <<  v1 << '\n';
          return true;
        } else {
          UnionFindUnionUtil(v1, v2, parent);
        }
      }
    }
  }
  return false;
}

int GraphAM::UnionFindRCFindUtil(const int v, std::vector<std::pair<int, int>>& subsets) const {
  if(subsets[v].first != v) {
    subsets[v].first = UnionFindRCFindUtil(subsets[v].first, subsets);
  }
  return subsets[v].first;
}

void GraphAM::UnionFindRCUnionUtil(const int v1, const int v2, std::vector<std::pair<int, int>>& subsets) const {
  int rv1 = UnionFindRCFindUtil(v1, subsets);
  int rv2 = UnionFindRCFindUtil(v2, subsets);

  if(subsets[rv1].second < subsets[rv2].second) {
    subsets[rv1].first =  rv2;
  } else if (subsets[rv1].second > subsets[rv2].second) {
    subsets[rv2].first = rv1;
  } else {
    subsets[rv2].first = rv1;
    subsets[rv1].second++;
  }
}

bool GraphAM::UnionFindRCDetectCycle() const {

  std::vector<std::pair<int, int>> subsets;
  subsets.reserve(v);
  for(int i=0;i<v;i++) {
    subsets.emplace_back(std::make_pair(i, 0));
  }
  for (int i = 0; i < v; i++) {
    for (int j = 0; j < v; j++) {
      if(g[i][j] != 0) {
        int v1 = UnionFindRCFindUtil(i, subsets);
        int v2 = UnionFindRCFindUtil(j, subsets);
        if(v1==v2) {
          std::cout << __FUNCTION__ << " | " <<  " Cycle detected" << '\n';
          std::cout << __FUNCTION__ << " | " <<  " Edge vertices " << i << ' ' << j << '\n';
          std::cout << __FUNCTION__ << " | " <<  " Common vertex " <<  v1 << '\n';
          return true;
        } else {
          UnionFindRCUnionUtil(v1, v2, subsets);
        }
      }
    }
  }
  return false;
}

void GraphAM::DFSUtil(const int source, std::vector<bool>& visited) const {

  visited[source] = true;
  for(int i=0; i < v; i++) {
    if(g[source][i] !=0  && !visited[i]) {
      DFSUtil(i, visited);
    }
  }
}

std::tuple<bool, int> GraphAM::FindMotherVertex() const {

  if(v == 0) return {false, -1};
  std::vector<bool> visited(v, false);
  int mv = -1;
  for(int i=0; i<v; i++) {
    if(!visited[i]) {
      DFSUtil(i, visited);
      mv = i;
    }
  }

  std::fill(visited.begin(), visited.end(), false);
  DFSUtil(mv, visited);
  for(const auto& ele : visited) {
    if(!ele) return {false, -1};
  }
  std::cout << __FUNCTION__ << " | " <<  " Mother vetex found" << '\n';
  std::cout << __FUNCTION__ << " | " <<  " Mother vertex is "<< mv << '\n';
  return {true, mv};
}

void GraphAM::allPathsBetweenUtil(const int source, int sink, std::vector<int>& path,
  std::vector<std::vector<int>>& paths, std::vector<bool>& visited) const {

  visited[sink] = true;
  path.push_back(sink);

  if(sink == source) {
    paths.push_back(path);
    path.pop_back();
    visited[sink] = false;
    return;
  }


  for(int i = 0; i < v ; i++) {
    if(g[sink][i] != 0 && !visited[i]) {
      allPathsBetweenUtil(source, i, path, paths, visited);
    }
  }

  path.pop_back();
  visited[sink] = false;
  return;
}

std::tuple<bool, std::vector<std::vector<int>>> GraphAM::allPathsBetween(const int source, int sink) const {

  if(v==0) return {false, std::vector<std::vector<int>>()};

  std::vector<int> path;
  std::vector<std::vector<int>> paths;
  std::vector<bool> visited(v, false);

  allPathsBetweenUtil(source, sink, path, paths, visited);
  if(paths.empty()) {
    std::cout << __FUNCTION__ << " | " <<  " No path found" << '\n';
    return {false, paths};
  } else {
    std::cout << __FUNCTION__ << " | " <<  ' ' << paths.size() << " paths found" << '\n';
    for(const auto& path : paths) {
      std::cout << __FUNCTION__ << " | " <<  " Path: ";
      for(const auto& vert : path) {
        std::cout << vert << ' ';
      }
      std::cout << '\n';
    }
    return {true, paths};
  }
}

bool GraphAM::colourGraphUtil(const int vert, const int n_c, std::vector<int>& colour) const {
  bool safe_to_colour = true;

  for(int c = 1; c <= n_c; c++) {
    colour[vert] = c;
    safe_to_colour = true;

    // Check if this is an acceptable colouring
    for(int i=0; i < v; i++) {
      if (colour[i] == c && (g[vert][i]!=0 || g[i][vert]!=0)) {
        safe_to_colour = false;
        break;
      }
    }
    if(!safe_to_colour) continue;

    // The recursive part of the algorithm
    for(int i = 0; i < v; i++) {
      if (colour[i] == 0 && (g[vert][i]!=0 || g[i][vert]!=0)) {  // can be directed graph
        safe_to_colour = colourGraphUtil(i, n_c, colour);
      }
      if(!safe_to_colour) break;
    }
    if(safe_to_colour) break;
  }

  if(!safe_to_colour) colour[vert] = 0;
  return safe_to_colour;
}

std::tuple<bool, std::vector<int>> GraphAM::colourGraph(const int n_c) const {

  if(v == 0 && n_c > 0 ) return {true, std::vector<int>()};
  if(v != 0 && n_c == 0 ) return {false, std::vector<int>()};

  for(int i = 0; i < v; i++) {
    if(g[i][i] != 0) {
      std::cout << __FUNCTION__ << " | " <<  " Vertex " << i << " has an edge to itelf " << '\n';
      std::cout << __FUNCTION__ << " | " <<  " Unable to colour graph using " << n_c << " colours " << '\n';
      return {false, std::vector<int>()};
    }
  }

  std::vector<int> colour(v, 0);

  for(int i = 0; i < v; i++) {
    if(colour[i] == 0) {
      bool coloured = colourGraphUtil(i, n_c, colour);
      if(!coloured) {
        std::cout << __FUNCTION__ << " | " <<  " Unable to colour graph using " << n_c << " colours " << '\n';
        return {false, std::vector<int>()};
      }
    }
  }
  std::cout << __FUNCTION__ << " | " <<  " Coloured the given graph using " << n_c << " colours " << '\n';
  return {true, colour};
}

// Cannot handle negative graph weights
std::tuple<bool, double, std::vector<int>> GraphAM::Dijkstra(const int source, const int sink) const {

  struct PQSortSecondPair { //
    bool operator()(const std::pair<int, double> p1, const std::pair<int, double> p2) {
      return p1.second > p2. second;
    }
  };


  std::vector<std::pair<int, double>> p_md(v, std::make_pair(-1, std::numeric_limits<double>::max())); // parent, min distance
  std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, PQSortSecondPair> pq; // <point id, cost to get there>
  pq.push({source, 0});
  bool found_path = false;
  std::pair<int, double> current = pq.top();
  while(!pq.empty()) {
    current = pq.top();
    pq.pop();
    if(current.first == sink) {
      found_path = true;
      break;
    }
    for(int i = 0; i < v; i++) {
      if((g[current.first][i]!=0) &&
      (current.second + g[current.first][i] < p_md[i].second)) {
        pq.push(std::make_pair(i, current.second + g[current.first][i]));
        p_md[i] = std::make_pair(current.first, current.second + g[current.first][i]);
      }
    }
  }

  if(found_path == false) {
    std::cout << __FUNCTION__ << "Path not found between " << source << ' ' << sink << '\n';
    return {false, -1, std::vector<int>()};
  } else {
    int cv = sink;
    double path_cost = 0;
    std::vector<int> path;
    while(cv!=source) {
      path.push_back(cv);
      cv = p_md[cv].first;
      path_cost += p_md[cv].second;
    }
    path.push_back(source);
    std::reverse(path.begin(), path.end());
    std::cout << __FUNCTION__ << "Path found between " << source << ' ' << sink << '\n';
    return {found_path, path_cost, path};
  }
}

void GraphAM::PrintAllDijkstraPathsFound(const int source, const std::vector<std::pair<int, double>>& p_md) const {
  for(int i = 0; i < p_md.size(); i++) {
    if(p_md[i].first != -1) {
      std::vector<int> path;
      path.push_back(i);
      std::cout << __FUNCTION__ <<  " | " << "Path: ";
      int cv = p_md[i].first;
      while(cv != source) {
        path.push_back(cv);
        cv = p_md[cv].first;
      }
      path.push_back(source);
      std::reverse(path.begin(), path.end());
      for(const auto& ele : path) {
        std::cout << ele << ' ';
      }
      std::cout << '\n';
      std::cout << __FUNCTION__ <<  " | " << "Cost: " << p_md[i].second << '\n';
    }
  }
}

std::tuple<bool, std::vector<std::pair<int, double>>> GraphAM::Dijkstra(const int source) const {

  struct PQSortSecondPair { //
    bool operator()(const std::pair<int, double> p1, const std::pair<int, double> p2) {
      return p1.second > p2. second;
    }
  };


  std::vector<std::pair<int, double>> p_md(v, std::make_pair(-1, std::numeric_limits<double>::max())); // parent, min distance
  std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, PQSortSecondPair> pq; // <point id, cost to get there>
  pq.push({source, 0});
  std::pair<int, double> current = pq.top();
  while(!pq.empty()) {
    current = pq.top();
    pq.pop();
    // Finding for all points so no early stop condition
    for(int i = 0; i < v; i++) {
      if((g[current.first][i]!=0) &&
      (current.second + g[current.first][i] < p_md[i].second)) {
        pq.push(std::make_pair(i, current.second + g[current.first][i]));
        p_md[i] = std::make_pair(current.first, current.second + g[current.first][i]);
      }
    }
  }
  // Check if any paths found
  bool found_any = false;
  for(const auto& ele : p_md) {
    if(ele.first!=-1) {
      found_any = true;
      break;
    }
  }
  PrintAllDijkstraPathsFound(source, p_md);
  return {found_any, p_md};
}

// visit
void GraphAM::DFSUtilWithFinishTime(int source, std::vector<bool>& visited, std::stack<int>& visit_order) const {
  visited[source] = true;
  for(int i = 0; i < v; i++) {
    if(!visited[i] && g[source][i]!=0) {
      DFSUtilWithFinishTime(i, visited, visit_order);
    }
  }
  visit_order.push(source);
}

std::vector<std::vector<int>> GraphAM::KosarajuAlgorithmUtil() {
  std::vector<bool> visited(v, false);
  std::stack<int> visit_order;
  for(int i = 0; i < v; i++) {
    if (!visited[i]) {
      DFSUtilWithFinishTime(i, visited, visit_order);
    }
  }

  for(int i=0; i< v; i++) {
    for(int j = i+1; j < v; j++) {
      std::swap(g[i][j], g[j][i]);
    }
  }

  std::fill(visited.begin(), visited.end(), false);

  std::vector<std::vector<int>> components;
  while(!visit_order.empty()) {
    int vert = visit_order.top();
    if(!visited[vert]) {
      std::stack<int> component_stack;
      std::vector<int> component_vector;
      DFSUtilWithFinishTime(vert, visited, component_stack);
      while(!component_stack.empty()) {
        component_vector.push_back(component_stack.top());
        component_stack.pop();
      }
      components.push_back(component_vector);
    }
    visit_order.pop();
  }
  return components;
}

// Needs to modify graph
std::vector<std::vector<int>> GraphAM::KosarajuAlgorithm() const {
  GraphAM g_am(g);
  return g_am.KosarajuAlgorithmUtil();
}

bool GraphAM::StronglyConnectedKosarajuUtil() {
  std::vector<bool> visited(v, false);
  std::stack<int> visit_order;
  DFSUtil(0, visited);

  for(const auto& ele : visited) {
    if(!ele) return false;
  }

  for(int i=0; i < v; i++) {
    for(int j = i+1; j< v; j++) {
      std::swap(g[i][j], g[j][i]);
    }
  }

  std::fill(visited.begin(), visited.end(), false);
  DFSUtil(0, visited);

  for(const auto& ele : visited) {
    if(!ele) return false;
  }

  return true;
}

bool GraphAM::StronglyConnectedKosaraju() const {
  GraphAM g_am(g);
  return g_am.StronglyConnectedKosarajuUtil();
}

std::vector<Edge> GraphAM::KruskalsAlgorithm() const {

  struct EdgeCompare {
    bool operator() (const Edge e1, const Edge e2) const {
      return e1.w > e2.w;
    }
  };

  std::priority_queue<Edge, std::vector<Edge>, EdgeCompare> pq;
  for(int i=-0; i<v;i++) {
    for(int j = 0; j<v; j++){
      if(g[i][j]!=0) {
        pq.push(Edge(i, j, g[i][j]));
      }
    }
  }

  std::vector<Edge> mst;
  // std::vector<int> parents(v, -1);
  std::vector<std::pair<int,int>> subsets;
  subsets.reserve(v);
  for(int i=0; i< v; i++) {
    subsets.emplace_back(i,0);
  }

  while(mst.size() < v-1 && !pq.empty()) {
    auto e = pq.top();
    pq.pop();
    if(UnionFindRCFindUtil(e.u, subsets)!=UnionFindRCFindUtil(e.v, subsets)) {
      mst.push_back(e);
      UnionFindRCUnionUtil(e.u,e.v, subsets);
    }
  }
  return mst;
}

std::vector<std::vector<double>> GraphAM::FloydWarshall() const {
  std::vector<std::vector<double>> fwg = g;
  for(int k = 0; k < v; k++) {
    for(int i = 0; i < v; i++) {
      for(int j = 0; j < v; j++) {
        if(i!=j && (fwg[i][j] == 0 || fwg[i][k] + fwg[k][j] < fwg[i][j] )) {  // 0 implies no edge
          fwg[i][j] = fwg[i][k] + fwg[k][j];
        }
      }
    }
  }
  return fwg;
}

// TODO: To test
bool GraphAM::BellmanFord(const int source) const {
  std::vector<int> parents(v, -1);
  std::vector<double> distance(v, std::numeric_limits<double>::max());
  distance[source] = 0;
  // NOTE: As of now a weight of 0 implies no edge exists between 2 points
  // Converting the graph to std::numeric_limits<double>::max() to represent
  // no edge and accepting 0 as a possible edge weight for this algorithm
  auto g_with_negative = g;
  for(int j = 0; j < v; j++) {
    for(int k = 0; k < v; k++) {
      g_with_negative[j][k] = g[j][k];
      if(!g_with_negative[j][k]) {
        g_with_negative[j][k] = std::numeric_limits<double>::max();
      }
    }
  }

  for(int i = 0; i < v-1; i++) {
    for(int j = 0; j < v; j++) {
      for(int k = 0; k < v; k++) {
        if(distance[k] > distance[j] + g_with_negative[j][k]) {
          distance[k] =  distance[j] + g_with_negative[j][k];
          parents[k] = j;
        }
      }
    }
  }

  for(int j = 0; j < v; j++) {
    for(int k = 0; k < v; k++) {
      if(distance[k] > distance[j] + g_with_negative[j][k]) {
        return true;
      }
    }
  }
  return false;
}

} // namespace graphAM


// int main () {
//
//   // Common graph for most of the tests below
//   graphAM::GraphAM g_am;
//
//   // Test random distribution
//   std::binomial_distribution<> dist(1,0.2);
//   graphAM::GraphAM g1(11, dist);
//   std::cout << "Test random distribution " << '\n';
//   g1.PrintGraph();
//
//   // Random with different weights
//   graphAM::GraphAM g2(11, true, true);
//   std::cout << "Test graph with different_weights " << '\n';
//   g2.PrintGraph();
//
//   // Random with same weights
//   graphAM::GraphAM g3(11, true, false);
//   std::cout << "Test graph with same weights " << '\n';
//   g3.PrintGraph();
//
//   // Empty graph
//   graphAM::GraphAM g4(11, false);
//   std::cout << "Test empty graph " << '\n';
//   g4.PrintGraph();
//
//   // No argumenmts specified
//   graphAM::GraphAM g5;
//   std::cout << "Test no arguments specified " << '\n';
//   g5.PrintGraph();
//
//   // Test for BFS/DFS/Prim
//   std::vector<std::vector<double>> g = {
//     { 1, 1, 1, 0, 0 },
//     { 0, 1, 1, 0, 0 },
//     { 0, 0, 1, 1, 0 },
//     { 0, 0, 0, 1, 1 },
//     { 0, 0, 0, 0, 1 },
//   };
//
//   g_am = graphAM::GraphAM(g);
//   auto [path_found_bfs, path_bfs] = g_am.BFS(0, 4);
//   if(path_found_bfs){
//     std::cout << "Path: " << '\n';
//     for(const auto& ele : path_bfs) std::cout << ele << "   ";
//     std::cout << '\n';
//   } else {
//     std::cout << "No path found" << '\n';
//   }
//
//   auto [path_found_dfs, path_dfs] = g_am.DFS(0, 4);
//   if(path_found_dfs){
//     std::cout << "Path: " << '\n';
//     for(const auto& ele : path_dfs) std::cout << ele << "   ";
//     std::cout << '\n';
//   } else {
//     std::cout << "No path found" << '\n';
//   }
//
//   // Test for Prim
//   auto found_prim = g_am.Prim();
//
//   g = {
//     {0, 16, 13, 0, 0, 0},
//     {0, 0, 10, 12, 0, 0},
//     {0, 4, 0, 0, 14, 0},
//     {0, 0, 9, 0, 0, 20},
//     {0, 0, 0, 7, 0, 4},
//     {0, 0, 0, 0, 0, 0}
//   };
//
//   g_am = graphAM::GraphAM(g);
//   double flow = g_am.fordFulkerson(0, 5);
//   std::cout << "Flow from source to sink: " << flow << '\n';
//
//   // Test for mother vertex
//   g = {
//     {0, 1, 1, 0, 0, 0, 0},
//     {0, 0, 0, 1, 0, 0, 0},
//     {0, 0, 0, 0, 0, 0, 0},
//     {0, 0, 0, 0, 0, 0, 0},
//     {0, 1, 0, 0, 0, 0, 0},
//     {0, 0, 1, 0, 0, 0, 1},
//     {1, 0, 0, 0, 1, 0, 0}
//   };
//
//   g_am = graphAM::GraphAM(g);
//   auto [mother_vertex_found, mv] = g_am.FindMotherVertex();
//
//   // Test for UnionFind and UnionFindRC
//   g = {
//     {0, 1, 1, 0, 0, 0},
//     {0, 0, 0, 1, 0, 0},
//     {0, 0, 0, 0, 1, 0},
//     {0, 0, 0, 0, 0, 1},
//     {0, 0, 0, 0, 0, 1},
//     {0, 0, 0, 0, 0, 0}
//   };
//
//   g_am = graphAM::GraphAM(g);
//   bool cycle_found = g_am.UnionFindRCDetectCycle();
//   std::cout << "Cycle found in graph: " << cycle_found << '\n';
//
//   // Test for allPathsBetween
//   g = {
//     {0, 1, 1, 0, 0, 0},
//     {1, 0, 1, 0, 0, 1},
//     {1, 1, 0, 1, 1, 1},
//     {0, 0, 1, 0, 1, 1},
//     {0, 0, 1, 1, 0, 1},
//     {0, 1, 1, 1, 1, 0}
//   };
//
//   g_am = graphAM::GraphAM(g);
//   auto [paths_found, paths] = g_am.allPathsBetween(0, 5);
//
//   // Test for colourGraph
//   std::vector<std::vector<double>> g = {
//     {0, 1, 1, 0, 0, 0},
//     {1, 0, 1, 0, 0, 1},
//     {1, 1, 0, 1, 1, 1},
//     {0, 0, 1, 0, 1, 1},
//     {0, 0, 1, 1, 0, 1},
//     {0, 1, 1, 1, 1, 0}
//   };
//
//   g_am = graphAM::GraphAM(g);
//   auto [coloured, colours] = g_am.colourGraph(4);
//
// g = {
//    {0, 1, 1, 0, 0, 0},
//    {1, 0, 1, 0, 0, 0},
//    {1, 1, 0, 1, 1, 0},
//    {0, 0, 1, 0, 1, 0},
//    {0, 0, 1, 1, 0, 1},
//    {0, 1, 1, 1, 1, 0}
//   };
//
//   g_am = graphAM::GraphAM(g);
//   // const auto[path_found, path_cost, path] = g_am.Dijkstra(0, 5);
//   const auto[path_found, p_md] = g_am.Dijkstra(0);
//
//   g = {
//    {0, 0, 1, 1, 0},
//    {1, 0, 0, 0, 0},
//    {0, 1, 0, 0, 0},
//    {0, 0, 0, 0, 1},
//    {0, 0, 0, 0, 0},
//   };
//   g_am = graphAM::GraphAM(g);
//   const auto components = g_am.KosarajuAlgorithm();
//   for(auto& row : components) {
//     std::cout << "Connected component: ";
//     for(auto & element : row) {
//       std::cout << element;
//     }
//     std::cout << '\n';
//   }
//
//   g = {
//    {0, 0, 1, 1, 0},
//    {1, 0, 0, 0, 0},
//    {1, 1, 0, 0, 0},
//    {1, 0, 0, 0, 1},
//    {0, 0, 0, 1, 0},
//   };
//   g_am = graphAM::GraphAM(g);
//   const auto strongly_connected = g_am.StronglyConnectedKosaraju();
//   if(strongly_connected) {
//     std::cout << "The graph is strongly connected" << '\n';
//   } else {
//     std::cout << "The graph is not strongly connected" << '\n';
//   }
//
//   return 0;
// }
