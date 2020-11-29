#include <algorithm>
#include <iostream>
#include <queue>
#include <limits>
#include <random>
#include <stack>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include <chrono>
#include <thread>

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
        std::cout << "Edge from " << i << " to " << j <<
        " not symetric/undirected. Prim's algorithm  will fail. " << '\n';
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
      if (i!=source) edge_pq.push({source, i, g[source][i]});
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
          edge_pq.push({e.v, i, g[e.v][i]});
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
    path_cost = p_md[cv].second;
    while(cv!=source) {
      path.push_back(cv);
      cv = p_md[cv].first;
    }
    path.push_back(source);
    std::reverse(path.begin(), path.end());
    std::cout << __FUNCTION__ << "Path found between " << source << ' ' << sink << '\n';
    std::cout << __FUNCTION__ << "Path cost: " << path_cost << '\n';
    std::cout << __FUNCTION__ << "Path: ";
    // for(const auto & p : path) {
    //   std::cout << p << ' ';
    // }
    // std::cout << '\n';
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
        pq.push({i, j, g[i][j]});
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

std::vector<int> GraphAM::HierholzersAlgorithmUtil() {
  // check degrees
  std::vector<int> indeg(v, 0), outdeg(v, 0);
  for(int i = 0 ; i < v; i++) {
    for(int j = 0; j < v; j++) {
      if (g[i][j] != 0) {
        ++indeg[j];
        ++outdeg[i];
      }
    }
  }

  // start vertex
  int first  = 0;
  while(first < v && !indeg[first]) {
    ++ first;
  }
  if(first == v) {
    std::cout << "The are no edges in this graph" << '\n';
    return std::vector<int>();
  }

  // check if at maximum only 2 vertices have indeg != outdeg and if so, connect them

  int v1 = -1, v2 = -1, diff = 0;
  for(int i = 0; i < v; i++) {
    if(indeg[i] - outdeg[i] != 0) {
      if(v1 == -1) {
        v1 = i;
      } else if (v2 == -1) {
        v2 = i;
      }
      else {
        std::cout << "There exist more than 2 vertices for which the in degree is not equal to the out degree" << '\n';
        return std::vector<int>();
      }
    }
  }

  if(v1 != -1) {
    if(indeg[v1] - outdeg[v1] > 0) {
      g[v1][v2] = 1;
    } else {
      g[v2][v1] = 1;
    }
  }

  std::stack<int> vertices; // vertices to expand
  std::vector<int> eulerian_path;
  vertices.push(first);

  while(!vertices.empty()) {
    const int edge_start = vertices.top();
    int edge_end;
    for(edge_end = 0; edge_end < v; edge_end++) {
      if(g[edge_start][edge_end] != 0) {
        break;
      }
    }
    if(edge_end == v) {
      eulerian_path.push_back(edge_start);
      vertices.pop();
    } else {
      g[edge_start][edge_end] = 0;
      vertices.push(edge_end);
    }
  }

  // remove added edge if added
  if(v1 != -1) {
    const int s1 = eulerian_path.size() - 1;
    for(int i = 0; i < s1; i++) {
      if((eulerian_path[i] == v1 && eulerian_path[i+1] == v2) ||
         (eulerian_path[i] == v2 && eulerian_path[i+1] == v1)) {
        eulerian_path.pop_back();
        std::rotate(std::begin(eulerian_path), std::next(std::begin(eulerian_path), i), std::end(eulerian_path));
        break;
      }
    }
  }

  // check whether any edge remains
  // TODO: check when this is required
  for(const auto& row:g) {
    for(const auto& ele : row) {
      if(ele != 0) {
        std::cout << "There is an edge left" << '\n';
        return std::vector<int>();
      }
    }
  }

  std::reverse(eulerian_path.begin(), eulerian_path.end());
  return eulerian_path;
}

std::vector<int> GraphAM::HierholzersAlgorithm() const {
  GraphAM mod(g);
  return mod.HierholzersAlgorithmUtil();
}

void GraphAM::ArticulationPointsUtil(const int vert1, std::vector<bool>& visited, std::vector<int>& parent,
  std::vector<int>& tod, std::vector<int>& low, std::unordered_set<int>& articulation_points, int time) const {
    int children = 0;
    visited[vert1] = true;
    tod[vert1] = low[vert1] = ++time;
    for(int vert2 = 0; vert2 < v; ++vert2) {
      if(g[vert1][vert2] == 0) continue; // check if this should also check for parent
      if(!visited[vert2]) {
        ++children;
        parent[vert2] = vert1;
        ArticulationPointsUtil(vert2, visited, parent, tod, low, articulation_points, time);
        low[vert1] = std::min(low[vert1], low[vert2]);
        if(low[vert2] >= tod[vert1] && parent[vert1] != -1) {
          articulation_points.insert(vert1);
        }
      } else {
        low[vert1] = std::min(low[vert1], tod[vert2]);
      }
    }
    if(parent[vert1] == -1 && children > 1) {
      articulation_points.insert(vert1);
    }
}


std::unordered_set<int> GraphAM::ArticulationPoints() const {

  // While this algorithm will correctly find all teh articulation points as long as
  // for each edge between 2 vertices there exists an edge runnin in the opposite direction
  // regardless of weight, this algorithm returns when the graph is not symetric to emphasize
  // the fact that this algorithm is meant to be run n undirected graphs
  for(int i = 0; i < v; i++) {
    for (int j = i+1; j < v; j++) {
      if(g[i][j] != g[j][i]) {
        std::cout << "Edge from " << i << " to " << j <<
        " not symetric. The algorithm will fail. " << '\n';
        return std::unordered_set<int>();
      }
    }
  }

  std::vector<bool> visited(v, false);
  std::vector<int> parent(v, -1);
  std::vector<int> tod(v, std::numeric_limits<int>::max()); // time of discovery
  std::vector<int> low(v, std::numeric_limits<int>::max());
  std::unordered_set<int> articulation_points;
  int time = 0;

  for(int i = 0; i < v; ++i) {
    if(!visited[i]) {
      ArticulationPointsUtil(i, visited, parent, tod, low, articulation_points, time);
    }
  }
  return articulation_points;
}

bool GraphAM::TopologicalSortUtil(std::vector<int>& sorted,
  std::vector<bool>& visited,
  std::vector<bool>& this_cycle,
  const int vert) const {

  this_cycle[vert] = true;
  visited[vert] = true;
  for(int i = 0; i < v; i++) {
    if(g[vert][i] != 0) {
      if(this_cycle[i] && visited[i]) {
        return false;
      } else if(!visited[i]) {
        if(!TopologicalSortUtil(sorted, visited, this_cycle, i)) {
          return false;
        }
      }
    }
  }
  this_cycle[vert] = false;
  sorted.push_back(vert);
  return true;
}

std::tuple<bool, std::vector<int>> GraphAM::TopologicalSort() const {
  std::vector<int> sorted;
  std::vector<bool> visited(v, false);
  std::vector<bool> this_cycle(v, false);
  for(int i = 0; i < v; ++i) {
    if(!visited[i]) {
      if(!TopologicalSortUtil(sorted, visited, this_cycle, i)) {
        return {false, std::vector<int>()};
      }
    }
  }
  std::reverse(sorted.begin(), sorted.end());
  return {true, sorted};
}

void GraphAM::FindBridgesUtil(const int vert1, std::vector<bool>& visited, std::vector<int>& parent,
  std::vector<int>& tod, std::vector<int>& low, std::vector<std::pair<int, int>>& bridges, int time) const {
  visited[vert1] = true;
  low[vert1] = tod[vert1] = time++;
  for(int vert2 = 0; vert2 < v; ++vert2) {
    if(g[vert1][vert2] == 0 || parent[vert1] == vert2) continue;
    if(!visited[vert2]) {
      parent[vert2] = vert1;
      FindBridgesUtil(vert2, visited, parent, tod, low, bridges, time);
      low[vert1] = std::min(low[vert1], low[vert2]);
      if(low[vert2] > tod[vert1]) {
        bridges.emplace_back(vert1, vert2);
      }
    } else {
      low[vert1] = std::min(low[vert1], tod[vert2]);
    }
  }
}

std::vector<std::pair<int, int>> GraphAM::FindBridges() const {

  for(int i = 0; i < v; i++) {
    for (int j = i+1; j < v; j++) {
      if(g[i][j] != g[j][i]) {
        std::cout << "Edge from " << i << " to " << j <<
        " not symetric. The algorithm will fail. " << '\n';
        return std::vector<std::pair<int, int>> ();
      }
    }
  }

  int time = 0;
  std::vector<bool> visited(v, false);
  std::vector<int> parent(v, -1);
  std::vector<int> tod(v, std::numeric_limits<int>::max());
  std::vector<int> low(v, std::numeric_limits<int>::max());
  std::vector<std::pair<int, int>> bridges;

  for(int i = 0; i < v; ++i) {
    if(!visited[i]) {
      FindBridgesUtil(i, visited, parent, tod, low, bridges, time);
    }
  }

  return bridges;
}

bool GraphAM::IsBipartite() const {
  std::vector<int> side(v, -1);
  std::queue<int> q;

  for(int vert1 = 0; vert1 < v; vert1++) {
    if(side[vert1] == -1) {
      q.push(vert1);
      side[vert1] = 0;
      while(!q.empty()) {
        const int vert2 = q.front();
        q.pop();
        for (int vert3 = 0; vert3 < v; ++vert3) {
          if(g[vert2][vert3] == 0) continue;
          if(side[vert3] == -1) {
            side[vert3] = side[vert2] ^ 1; // xor
            q.push(vert3);
          } else {
            if(side[vert2] == side[vert3]) {
              std::cout << vert2 << ' ' << vert3 << " have the same colour"<< '\n';
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}

bool GraphAM::DivideIntoTwoCliques() const {
  auto g2 = g;
  for(auto& row : g2) {
    for(auto& ele : row) {
      if(ele != 0) {
        ele = 0;
      } else {
        ele = 1;
      }
    }
  }
  for(int i = 0; i < v; i++) {
    g2[i][i] = 0;
  }

  GraphAM gcomp(g2);
  return gcomp.IsBipartite();
}

std::tuple<bool, std::vector<int>> GraphAM::CreateLevelGraph(const int source, const int sink) const {
  std::vector<int> levels(v, -1);
  int level = 0;
  std::queue<int> q;
  q.push(source);
  while(!q.empty()) {
    const int size = q.size();
    for(int i = 0; i < size; ++i) {
      const int vert = q.front();
      q.pop();
      levels[vert] = level;
      for(int vert2 = 0; vert2 < v; ++vert2) {
        if(g[vert][vert2] != 0 && levels[vert2] == -1) {
          q.push(vert2);
        }
      }
      ++level;
    }
  }
  return {levels[sink] != -1, levels};
}

std::tuple<bool, std::vector<int>> GraphAM::CreateLevelGraph(const int source, const int sink, const std::vector<std::vector<double>>& g) const {
  std::vector<int> levels(v, -1);
  int level = 0;
  std::queue<int> q;
  q.push(source);
  while(!q.empty()) {
    const int size = q.size();
    for(int i = 0; i < size; ++i) {
      const int vert = q.front();
      q.pop();
      if(levels[vert] != -1) continue;
      levels[vert] = level;
      for(int vert2 = 0; vert2 < v; ++vert2) {
        if(g[vert][vert2] != 0 && levels[vert2] == -1) {
          q.push(vert2);
        }
      }
    }
    ++level;
  }
  return {levels[sink] != -1, levels};
}

double GraphAM::SendDinacFlowUtil(const int source, const int sink, double flow, const std::vector<int>& levels, const std::vector<std::vector<double>>& max_flow, std::vector<std::vector<double>>& rem_capacity) {
  if(source == sink) return flow;
  for(int vert2 = 0; vert2 < v; ++vert2) {
    if(max_flow[source][vert2] != 0 &&
      rem_capacity[source][vert2] != 0 &&
      levels[vert2] == levels[source] + 1) {
      const double temp_flow = SendDinacFlowUtil(vert2, sink, std::min(flow, rem_capacity[source][vert2]), levels, max_flow, rem_capacity);
      if(temp_flow > 0) {
        g[source][vert2] += temp_flow;
        g[vert2][source] -= temp_flow;
        rem_capacity[source][vert2] -= temp_flow;
        rem_capacity[vert2][source] += temp_flow;
        return temp_flow;
      }
    }
  }
  return 0;
}

double GraphAM::DinacsAlgorithmUtil(const int source, const int sink, const std::vector<std::vector<double>>& original) {
  double flow = std::numeric_limits<int>::max();
  double total_flow = 0;
  const auto max_flow = original;
  auto rem_capacity = original;
  GraphAM g_cap(g);
  auto ans = CreateLevelGraph(source, sink, rem_capacity);
  while(std::get<0>(ans)) {
    total_flow += SendDinacFlowUtil(source, sink, flow, std::get<1>(ans), max_flow, rem_capacity);
    ans = CreateLevelGraph(source, sink, rem_capacity);
  }
  return total_flow;
}

double GraphAM::DinacsAlgorithm(const int source, const int sink) const {
  GraphAM residual(std::vector<std::vector<double>>(v, std::vector<double>(v,0)));
  return residual.DinacsAlgorithmUtil(source, sink, g);
}

} // namespace graphAM
