#include <algorithm>
#include <iostream>
#include <queue>
#include <limits>
#include <random>
#include <stack>
#include <tuple>
#include <utility>
#include <vector>

#include "graphAM.hpp"

namespace graphAM {

Edge::Edge(int u, int v, double w) {
  this->u = u;
  this->v = v;
  this->w = w;
}

bool Edge::operator < (const Edge& other) const {
  return w < other.w;
}
bool Edge::operator <= (const Edge& other) const {
  return w <= other.w;
}
bool Edge::operator > (const Edge& other) const {
  return w > other.w;
}
bool Edge::operator >= (const Edge& other) const {
  return w >= other.w;
}

void PrintGraph(const GraphAM& g) {
  std::cout << "GraphAM: " << '\n';
  for(const auto& row : g) {
    for(const auto ele : row) {
      std::cout << ele << ' ';
    }
    std::cout << '\n';
  }
}

GraphAM createGraphAM(const int V, const bool random, const bool different_weights) {
  GraphAM g = std::vector<std::vector<double>>(V, std::vector<double>(V,0));
  if(random) {
    if(different_weights){
      std::uniform_int_distribution<> dist(0, V);
      return createGraphAMUtil(V, dist);
    }
    else {
      std::uniform_int_distribution<> dist(0, 1);
      return createGraphAMUtil(V, dist);
    }
  }
  else {
    return GraphAM(V, std::vector<double>(V, 0));
  }
}

std::tuple<bool, std::vector<int>> BFS (GraphAM& g, int source , int sink) {
  int v = g.size();
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


std::tuple<bool, std::vector<int>> DFS (GraphAM& g, int source , int sink) {
  int v = g.size();
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
bool Prim(GraphAM& g) {
  // TODO: Add check if connected

  int v = g.size();
  if( v == 0 ) return false;

  // This can be converted to a map of <point, parent>
  // BUt given that by the end of this algorithm it will contain v elements
  // anyway, sticking to vctor for now
  // Vector size is unchanging so no copy operations; and access is O(1)
  std::vector<int> parent(v, -1);

  int source = 0;
  parent[source] = source;

  std::priority_queue<Edge> edge_pq;

  // Insert current vertex
  for (int i = 0; i< v; i++) {
    if( g[source][i] != 0) {
      if (i!=source) edge_pq.push(Edge(source, i, g[source][i]));
    }
  }

  // Iterate until no edges left
  while(!edge_pq.empty()) {
    Edge e = edge_pq.top();
    edge_pq.pop();
    if(parent[e.v] == -1){
      parent[e.v] = e.u;
      for (int i = 0; i< v; i++) {
        if( g[e.v][i] != 0 && parent[i] == -1 ) {
          edge_pq.push(Edge(e.v, i, g[e.v][i]));
        }
      }
    }
  }

  std::cout << "Edges: " << '\n';
  for(int i = 0; i < source ; i++) {
    std::cout << parent[i]  << " --- " << i << std::endl;
  }
  // Do not print source --- source
  for(int i = source + 1; i < v ; i++) {
    std::cout << parent[i]  << " --- " << i << std::endl;
  }

  for(const auto& pid: parent) {
    if(pid ==-1 ) {
      std::cout << "Point " << pid << " not  included in MST " << '\n';
      return false;
    }
  }

  return true;
}

double fordFulkerson(GraphAM& g, int source, int sink) {
  GraphAM rg(g);
  double max_flow = 0;
  while(true){
    auto [found_path, path] = BFS(rg, source , sink);
    if(!found_path) break;
    double flow = std::numeric_limits<double>::max();
    for(int i = 1; i < path.size(); i++ ) {
      flow = std::min(flow, rg[path[i-1]][path[i]]);
    }
    for(int i = 1; i < path.size(); i++ ) {
      rg[path[i-1]][path[i]] -= flow;
      rg[path[i]][path[i-1]] += flow;
    }
    max_flow +=flow;
  }
  return max_flow;
}


int UnionFindFindUtil(const int v, const std::vector<int>& parent) {
  if(parent[v] == -1) return v;
  return UnionFindFindUtil(parent[v], parent);
}

void UnionFindUnionUtil(int v1, int v2, std::vector<int>& parent) {
  // The two calls below to find are not required as in the union find algorithm,
  // the parents are already found and the unionutil is only called when they
  // are not equal. Left here for completeness (and future use) and as it's only going to be a
  // single call, not a recursive one as their parents will be -1
  int rv1 = UnionFindFindUtil(v1, parent);
  int rv2 = UnionFindFindUtil(v2, parent);
  if (rv1 != rv2) parent[rv2] = rv1;
}

bool UnionFindDetectCycle(const GraphAM& g) {
  int v = g.size();
  std::vector<int> parent(v, -1);
  for (int i = 0; i < v ; i++) {
    for (int j = 0; j < v; j++) {
      if (g[i][j] != 0) {
        int v1 = UnionFindFindUtil(i, parent);
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


int UnionFindRCFindUtil(int v, std::vector<std::pair<int, int>>& subsets) {
  if(subsets[v].first != v) {
    subsets[v].first = UnionFindRCFindUtil(subsets[v].first, subsets);
  }
  return subsets[v].first;
}

void UnionFindRCUnionUtil(int v1, int v2, std::vector<std::pair<int, int>>& subsets) {
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

bool UnionFindRCDetectCycle (GraphAM& g) {
  int v = g.size();
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
  return true;
}

void DFSUtil (const GraphAM& g, int source, std::vector<bool>& visited) {
  int v = g.size();
  visited[source] = true;
  for(int i=0; i < v; i++) {
    if(g[source][i] !=0) {
      DFSUtil(g, i, visited);
    }
  }
}

std::tuple<bool, int> FindMotherVertex (const GraphAM& g) {
  int v = g.size();
  if(v == 0) return {false, -1};
  std::vector<bool> visited(v, false);
  int mv = -1;
  for(int i=0; i<v; i++) {
    if(!visited[i]) {
      DFSUtil(g, i, visited);
      mv = i;
    }
  }

  std::fill(visited.begin(), visited.end(), false);
  DFSUtil(g, mv, visited);
  for(const auto& ele : visited) {
    if(!ele) return {false, -1};
  }
  std::cout << __FUNCTION__ << " | " <<  " Mother vetex found" << '\n';
  std::cout << __FUNCTION__ << " | " <<  " Mother vertex is "<< mv << '\n';
  return {true, mv};
}

void allPathsBetweenUtil(const GraphAM& g, int source, int sink, std::vector<int>& path,
  std::vector<std::vector<int>>& paths, std::vector<bool>& visited) {

  visited[sink] = true;
  path.push_back(sink);

  if(sink == source) {
    paths.push_back(path);
    path.pop_back();
    visited[sink] = false;
    return;
  }

  const int v = g.size();
  for(int i = 0; i < v ; i++) {
    if(g[sink][i] != 0 && !visited[i]) {
      allPathsBetweenUtil(g, source, i, path, paths, visited);
    }
  }

  path.pop_back();
  visited[sink] = false;
  return;
}

std::tuple<bool, std::vector<std::vector<int>>> allPathsBetween(const GraphAM& g, int source, int sink) {
  const int v = g.size();
  if(v==0) return {false, std::vector<std::vector<int>>()};

  std::vector<int> path;
  std::vector<std::vector<int>> paths;
  std::vector<bool> visited(v, false);

  allPathsBetweenUtil(g, source, sink, path, paths, visited);
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

bool colourGraphUtil(const GraphAM& g, int vert, const int n_c, std::vector<int>& colour) {
  bool safe_to_colour = true;
  const int v = g.size();
  for(int c = 1; c <= n_c; c++) {
    colour[vert] = c;
    safe_to_colour = true;

    // More efficient to iterate over all the elements once just to make sure
    // that there is no clear colour repition before running recursive searches
    // and colour assignments
    for(int i=0; i < v; i++) {
      if (g[vert][i]!=0 && colour[i] == c) {
        safe_to_colour = false;
        break;
      }
    }
    if(!safe_to_colour) continue;

    // The recursive part of the algorithm
    for(int i=0; i < v; i++) {
      if(g[vert][i]!=0){
        if(colour[i] == 0) {
          safe_to_colour = colourGraphUtil(g, i, n_c, colour);
        } else if (colour[i] == c) {  // Check to make sure any new colouring does not conflict
          safe_to_colour = false;       // with the current vertex's colour
        }
        if(!safe_to_colour) break;
      }
    }

    if(safe_to_colour) break;
  }

  if(!safe_to_colour) colour[vert] = 0;
  return safe_to_colour;
}

std::tuple<bool, std::vector<int>> colourGraph(const GraphAM& g, const int n_c) {
  int v = g.size();
  if(v == 0 && n_c > 0 ) return {true, std::vector<int>()};
  if(v != 0 && n_c == 0 ) return {false, std::vector<int>()};

  for(int i = 0 ; i< v; i++) {
    if(g[i][i] != 0) {
      std::cout << __FUNCTION__ << " | " <<  " Vertex " << i << " has an edge to itelf " << '\n';
      std::cout << __FUNCTION__ << " | " <<  " Unable to colour graph using " << n_c << " colours " << '\n';
      return {false, std::vector<int>()};
    }
  }

  std::vector<int> colour(v, 0);
  for(auto vert : colour) {
    if(vert == 0) {
      bool coloured = colourGraphUtil(g, vert, n_c, colour);
      if(!coloured) {
        std::cout << __FUNCTION__ << " | " <<  " Unable to colour graph using " << n_c << " colours " << '\n';
        return {false, std::vector<int>()};
      }
    }
  }
  std::cout << __FUNCTION__ << " | " <<  " Coloured the given graph using " << n_c << " colours " << '\n';
  return {true, colour};
}


} // namespace graphAM


// int main () {
//
//   // Test random distribution
//   std::binomial_distribution<> dist(1,0.2);
//   graphAM::GraphAM g1 = graphAM::createGraphAM(11, dist);
//   std::cout << "Test random distribution " << '\n';
//   graphAM::PrintGraph(g1);
//
//   // Random with different weights
//   graphAM::GraphAM g2 = graphAM::createGraphAM(11, true, true);
//   std::cout << "Test graph with different_weights " << '\n';
//   graphAM::PrintGraph(g2);
//
//   // Random with same weights
//   graphAM::GraphAM g3 = graphAM::createGraphAM(11, true, false);
//   std::cout << "Test graph with same weights " << '\n';
//   graphAM::PrintGraph(g3);
//
//   // Empty graph
//   graphAM::GraphAM g4 = graphAM::createGraphAM(11, false);
//   std::cout << "Test empty graph " << '\n';
//   graphAM::PrintGraph(g4);
//
//   // No argumenmts specified
//   graphAM::GraphAM g5 = graphAM::createGraphAM();
//   std::cout << "Test no arguments specified " << '\n';
//   graphAM::PrintGraph(g5);
//
//   // Test for BFS/DFS/Prim
//   graphAM::GraphAM g = {
//     { 1, 1, 1, 0, 0 },
//     { 0, 1, 1, 0, 0 },
//     { 0, 0, 1, 1, 0 },
//     { 0, 0, 0, 1, 1 },
//     { 0, 0, 0, 0, 1 },
//   };
//
//   auto [path_found_bfs, path_bfs] = graphAM::BFS(g, 0, 4);
//   if(path_found_bfs){
//     std::cout << "Path: " << '\n';
//     for(const auto& ele : path_bfs) std::cout << ele << "   ";
//     std::cout << '\n';
//   } else {
//     std::cout << "No path found" << '\n';
//   }
//
//   auto [path_found_dfs, path_dfs] = graphAM::DFS(g, 0, 4);
//   if(path_found_dfs){
//     std::cout << "Path: " << '\n';
//     for(const auto& ele : path_dfs) std::cout << ele << "   ";
//     std::cout << '\n';
//   } else {
//     std::cout << "No path found" << '\n';
//   }
//
//   // Test for Prim
//   auto found_prim = graphAM::Prim(g);
//
// graphAM::GraphAM g = {
//   {0, 16, 13, 0, 0, 0},
//   {0, 0, 10, 12, 0, 0},
//   {0, 4, 0, 0, 14, 0},
//   {0, 0, 9, 0, 0, 20},
//   {0, 0, 0, 7, 0, 4},
//   {0, 0, 0, 0, 0, 0}
// };
//
// double flow = graphAM::fordFulkerson(g, 0, 5);
// std::cout << "Flow from source to sink: " << flow << '\n';
//
// return 0;
//   return 0;
// }

// Test for mother vertex
// graphAM::GraphAM g = {
//   {0, 1, 1, 0, 0, 0, 0},
//   {0, 0, 0, 1, 0, 0, 0},
//   {0, 0, 0, 0, 0, 0, 0},
//   {0, 0, 0, 0, 0, 0, 0},
//   {0, 1, 0, 0, 0, 0, 0},
//   {0, 0, 1, 0, 0, 0, 1},
//   {1, 0, 0, 0, 1, 0, 0}
// };
//
// auto [mother_vertex_found, mv] = graphAM::FindMotherVertex(g);
//
// Test for UnionFind and UnionFindRC
// graphAM::GraphAM g = {
//   {0, 1, 1, 0, 0, 0},
//   {0, 0, 0, 1, 0, 0},
//   {0, 0, 0, 0, 1, 0},
//   {0, 0, 0, 0, 0, 1},
//   {0, 0, 0, 0, 0, 1},
//   {0, 0, 0, 0, 0, 0}
// };
//
// bool cycle_found = graphAM::UnionFindRCDetectCycle(g);
// std::cout << "Cycle found in graph: " << cycle_found << '\n';
//
// Test for allPathsBetween
// graphAM::GraphAM g = {
//   {0, 1, 1, 0, 0, 0},
//   {1, 0, 1, 0, 0, 1},
//   {1, 1, 0, 1, 1, 1},
//   {0, 0, 1, 0, 1, 1},
//   {0, 0, 1, 1, 0, 1},
//   {0, 1, 1, 1, 1, 0}
// };
//
// auto [paths_found, paths] = graphAM::allPathsBetween(g, 0, 5);
// 
// Test for colourGraph
// graphAM::GraphAM g = {
//   {0, 1, 1, 0, 0, 0},
//   {1, 0, 1, 0, 0, 1},
//   {1, 1, 0, 1, 1, 1},
//   {0, 0, 1, 0, 1, 1},
//   {0, 0, 1, 1, 0, 1},
//   {0, 1, 1, 1, 1, 0}
// };
//
// auto [paths_found, paths] = graphAM::colourGraph(g, 4);
