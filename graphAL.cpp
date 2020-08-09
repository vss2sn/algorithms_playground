#include <iostream>
#include <queue>
#include <stack>
#include <utility>
#include <unordered_set>
#include <vector>

#include "graphAL.hpp"

namespace graphAL {

using Pair = std::pair<int, double>;

void GraphAL::PrintGraph() const {
  for(int i = 0; i < v; i++) {
    std::cout << "Vertex " << i << " ---> ";
    for(const auto& edge : g[i]) {
      std::cout << "(" << edge.first << ", " << edge.second << ") ";
    }
    std::cout << '\n';
  }
}

GraphAL::GraphAL(const std::vector<std::vector<Pair>>& g_v) {
  g = g_v;
  v = g.size();
}

std::tuple<bool, std::vector<int>> GraphAL::BFS(int source, int sink) const {

  if(source > v || sink > v)  return {false, std::vector<int>()};
  if(source == sink) return {true, std::vector<int>()};
  std::queue<int> q;
  std::vector<int> parent(v, -1);
  q.push(source);

  bool found = false;
  while(!q.empty()) {
    int vert = q.front();
    q.pop();
    for(const auto& p : g[vert]) {
      // the p.second!=0 check is rewuried for for fulkerson
      if(parent[p.first] == -1 && p.second != 0) {
        parent[p.first] = vert;
        if(p.first == sink) {
          found = true;
          break;
        }
        q.push(p.first);
      }
    }
  }
  if(found) {
    int vert_id = sink;
    std::vector<int> path;
    path.push_back(sink);
    while(vert_id!=source) {
      vert_id = parent[vert_id];
      path.push_back(vert_id);
    }
    std::reverse(path.begin(), path.end());
    return {true, path};
  } else {
    return {false, std::vector<int>()};
  }
}

std::tuple<bool, std::vector<int>> GraphAL::DFS(int source, int sink) const {
  if(source > v || sink > v) return {false, std::vector<int>()};
  if(source == sink) return {true, std::vector<int>()};

  std::stack<int> vert_stack;
  vert_stack.push(source);

  std::vector<int> parent(v, -1);
  bool found = false;

  while(!vert_stack.empty()) {
    int vert = vert_stack.top();
    vert_stack.pop();
    for(const auto& p : g[vert]) {
      if(parent[p.first] == -1) {
        parent[p.first] = vert;
        if(p.first == sink) {
          found = true;
          break;
        }
        vert_stack.push(p.first);
      }
    }
  }

  if(found) {
    int vert_id = sink;
    std::vector<int> path;
    path.push_back(sink);
    while(vert_id!=source) {
      vert_id = parent[vert_id];
      path.push_back(vert_id);
    }
    std::reverse(path.begin(), path.end());
    return {true, path};
  } else {
    return {false, std::vector<int>()};
  }
}

std::tuple<bool, std::vector<Edge>> GraphAL::Prim() const {
  if( v == 0 ) return {false, std::vector<Edge>()};

  // Check for duplicate edges between pairs of points
  for(int i = 0; i < v; i++) {
    std::unordered_set<int> end_points;
    for(const auto& [end_vert, weight] : g[i]) {
      const auto [iter, inserted] = end_points.insert(end_vert);
      if(!inserted) {
        std::cout << "Multiple edges from " << i << " to " << end_vert << '\n';
        std::cout << "Prim's algorithm will fail" << '\n';
        return  {false, std::vector<Edge>()};
      }
    }
  }

  std::vector<int> parent(v, -1);
  struct EdgeCompare {
    bool operator () (const Edge& e1, const Edge& e2) const {
      return e1.w > e2.w;
    }
  };

  std::priority_queue<Edge, std::vector<Edge>, EdgeCompare> pq;
  int source  = 0;
  parent[source] = source;

  for(const auto& [end_vert, weight] : g[source]) {
    pq.push(Edge(source, end_vert, weight));
  }

  std::vector<Edge> mst;
  while(!pq.empty() && mst.size() < v) {
    Edge current_edge = pq.top();
    pq.pop();
    if(parent[current_edge.v] == -1 ) {
      parent[current_edge.v] = current_edge.u;
      mst.push_back(current_edge);
      for(const auto& [end_vert, weight] : g[current_edge.v]) {
        pq.push(Edge(current_edge.v, end_vert, weight));
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

int GraphAL::UnionFindFindUtil(const int v, const std::vector<int>& parent) const {
  if(parent[v] == -1) return v;
  return UnionFindFindUtil(parent[v], parent);
}

void GraphAL::UnionFindUnionUtil(int v1, int v2, std::vector<int>& parent) const {
  int p1 = UnionFindFindUtil(v1, parent);
  int p2 = UnionFindFindUtil(v2, parent);
  if(p1 != p2) parent[p2] = p1;
}

bool GraphAL::UnionFindDetectCycle() const {
  std::vector<int> parent(v, -1);

  for(int start_vert = 0; start_vert < v ; start_vert++) {
    int v1 = UnionFindFindUtil(start_vert, parent);
    for(const auto& [end_vert, weight] : g[start_vert]) {
      if(start_vert != end_vert) {
        int v2 = UnionFindFindUtil(end_vert, parent);
        if(v1 == v2) {
          std::cout << __FUNCTION__ << " | " <<  " Cycle detected" << '\n';
          std::cout << __FUNCTION__ << " | " <<  " Edge vertices " << start_vert << ' ' << end_vert << '\n';
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

std::vector<Edge> GraphAL::KruskalsAlgorithm() const {

  std::vector<int> parent(v, -1);

  struct EdgeCompare {
    bool operator () (const Edge& e1, const Edge& e2) {
      return e1.w > e2.w;
    }
  };

  std::priority_queue<Edge, std::vector<Edge>, EdgeCompare> edge_pq;
  for(int start_vert = 0; start_vert < v; start_vert++) {
    for(const auto& [end_vert, weight] : g[start_vert]) {
      edge_pq.push(Edge(start_vert, end_vert, weight));
    }
  }

  std::vector<Edge> mst;
  while(mst.size() < v - 1 && !edge_pq.empty()) {
    Edge current_edge =  edge_pq.top();
    edge_pq.pop();
    if (UnionFindFindUtil(current_edge.u, parent) != UnionFindFindUtil(current_edge.v, parent)) {
      mst.push_back(current_edge);
      UnionFindUnionUtil(current_edge.u, current_edge.v, parent);
    }
  }
  return mst;
}

int GraphAL::UnionFindRCFindUtil(int v, std::vector<std::pair<int, int>>& subsets) const{
  if(v != subsets[v].first) {
    subsets[v].first = UnionFindRCFindUtil(subsets[v].first, subsets);
  }
  return subsets[v].first;
}

void GraphAL::UnionFindRCUnionUtil(int v1, int v2, std::vector<std::pair<int, int>>& subsets) const {
  int p1 = UnionFindRCFindUtil(v1, subsets);
  int p2 = UnionFindRCFindUtil(v2, subsets);

  if(subsets[p1].second > subsets[p2].second) {
    subsets[p2].first = p1;
  } else if(subsets[p1].second < subsets[p1].second) {
    subsets[p1].first = p2;
  } else {
    subsets[p2].first = p1;
    subsets[p1].second++;
  }
}

bool GraphAL::UnionFindRCDetectCycle() const {
  std::vector<std::pair<int, int>> subsets(v);
  for(int i=0; i < v; i++) {
    subsets.emplace_back(std::make_pair(i, 0));
  }

  for(int start_vert = 0; start_vert < v; start_vert) {
    int v1 = UnionFindRCFindUtil(start_vert, subsets);
    for(const auto& [end_vert, weights] : g[start_vert]) {
      int v2 = UnionFindRCFindUtil(end_vert, subsets);
      if(v1 == v2) {
        std::cout << __FUNCTION__ << " | " <<  " Cycle detected" << '\n';
        std::cout << __FUNCTION__ << " | " <<  " Edge vertices " << start_vert << ' ' << end_vert << '\n';
        std::cout << __FUNCTION__ << " | " <<  " Common vertex " <<  v1 << '\n';
        return true;
      } else {
        UnionFindRCUnionUtil(v1, v2, subsets);
      }
    }
  }
  return false;
}

double GraphAL::fordFulkersonRGUtil(int source, int sink) {
  double max_flow = 0;
  while(true) {
    auto [found_path, path] = BFS(source, sink);
    if(!found_path) {
      break;
    }
    double min_flow = std::numeric_limits<double>::max();
    for(int i = 0; i < path.size() - 1; i++) {
      for(const auto& [end_vert, weight] : g[path[i]]) {
        if(end_vert == path[i+1]) {
          min_flow = std::min(min_flow, weight);
        }
      }
    }
    for(int i = 0; i < path.size() - 1; i++) {
      for(auto& [end_vert, weight] : g[path[i]]) {
        if(end_vert == path[i+1]) {
          weight -= min_flow;
        }
      }
      for(auto& [end_vert, weight] : g[path[i+1]]) {
        if(end_vert == path[i]) {
          weight += min_flow;
        }
      }
    }
    max_flow += min_flow;
  }
  return max_flow;
}

double GraphAL::fordFulkerson(int source, int sink) const {
  GraphAL rg(*this);
  return rg.fordFulkersonRGUtil(source, sink);
}

void GraphAL::DFSUtil(const int v, std::vector<bool>& visited) const {
  visited[v] = true;
  for(const auto& [end_vert, weight] : g[v]) {
    if(!visited[end_vert] && weight != 0) DFSUtil(end_vert, visited);
  }
}

std::tuple<bool, int> GraphAL::FindMotherVertex () const {
  if(v == 0) return {false, -1};
  std::vector<bool> visited(v, false);
  int mv = -1;
  for(int i = 0; i < v; i++) {
    if(!visited[i]) {
      DFSUtil(i, visited);
      mv = i;
    }
  }
  std::fill(visited.begin(), visited.end(), false);
  DFSUtil(mv, visited);
  for(const auto& ele : visited) {
    if(!ele) return {false,-1};
  }
  std::cout << __FUNCTION__ << " | " <<  " Mother vetex found" << '\n';
  std::cout << __FUNCTION__ << " | " <<  " Mother vertex is "<< mv << '\n';
  return {true, mv};
}

std::tuple<bool, double, std::vector<int>> GraphAL::Dijkstra(const int source, const int sink) const {
  // point, min distance
  std::vector<std::pair<int, double>> p_md(v, std::make_pair(-1, std::numeric_limits<double>::max())); // parent, min distance
  struct PairSortSecond {
    bool operator () (const std::pair<int, double>& p1, const std::pair<int, double>& p2) {
      return p1.second > p2.second;
    }
  };

  std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, PairSortSecond> pq;

  p_md[source].second = 0;
  pq.push(std::make_pair(source,0));

  bool found_path = false;
  while(!pq.empty()) {
    const auto [current_vert, current_cost] = pq.top();
    pq.pop();
    if(current_vert == sink) {
      found_path = true;
      break;
    }
    if(p_md[current_vert].first == -1) {
      for(const auto& [end_vert, weight] : g[current_vert]) {
        if(double cost = current_cost + weight; cost < p_md[end_vert].second) {
          pq.push(std::make_pair(end_vert, cost));
          p_md[end_vert] = std::make_pair(current_vert, cost);
        }
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

void GraphAL::PrintAllDijkstraPathsFound(const int source, const std::vector<std::pair<int, double>>& p_md) const {
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

std::tuple<bool, std::vector<std::pair<int, double>>> GraphAL::Dijkstra(const int source) const {
  // point, min distance
  std::vector<std::pair<int, double>> p_md(v, std::make_pair(-1, std::numeric_limits<double>::max())); // parent, min distance
  struct PairSortSecond {
    bool operator () (const std::pair<int, double>& p1, const std::pair<int, double>& p2) {
      return p1.second > p2.second;
    }
  };

  std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, PairSortSecond> pq;

  p_md[source].second = 0;
  pq.push(std::make_pair(source,0));

  while(!pq.empty()) {
    const auto [current_vert, current_cost] = pq.top();
    pq.pop();
    if(p_md[current_vert].first == -1) {
      for(const auto& [end_vert, weight] : g[current_vert]) {
        if(double cost = current_cost + weight; cost < p_md[end_vert].second) {
          pq.push(std::make_pair(end_vert, cost));
          p_md[end_vert] = std::make_pair(current_vert, cost);
        }
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

bool GraphAL::colourGraphUtil(const std::vector<std::vector<std::pair<int, double>>>& g_undirected, const int vert, const int n_c, std::vector<int>& colour) const {
  std::cout << __FUNCTION__ << __LINE__ << '\n';
  std::cout << "Vert: " << vert << '\n';
  bool safe_to_colour = true;
  for(const auto& c: colour ) {
    std::cout << c << " ";
  }
  std::cout << '\n';
  for(int c = 1; c <= n_c; c++) {
    colour[vert] = c;
    safe_to_colour = true;
    for(const auto& [end_vert, weight] : g_undirected[vert]) {
      if(colour[end_vert] == c) {
        safe_to_colour = false;
        break;
      }
    }
    if(!safe_to_colour) {
      continue;
    }
    for(const auto& [end_vert, weight] : g_undirected[vert]) {
      if(colour[end_vert] == 0) {
        std::cout << __FUNCTION__ << __LINE__ << '\n';
        safe_to_colour = colourGraphUtil(g_undirected, end_vert, n_c, colour);
        if(!safe_to_colour) std::cout << "returned false" << '\n';
      }
      if(!safe_to_colour) break;
    }
    if(safe_to_colour) break;
  }
  if(!safe_to_colour) {
    colour[vert] = 0;
  }
  return safe_to_colour;
}

std::tuple<bool, std::vector<int>> GraphAL::colourGraph(int n_c) const {
  if(v == 0 && n_c > 0 ) return {true, std::vector<int>()};
  if(v != 0 && n_c == 0 ) return {false, std::vector<int>()};

  // Convert to undirected graph
  // TODO (vss): Make this more efficient and add as util
  std::vector<std::vector<std::pair<int, double>>> g_undirected = g;
  for(int i=0; i < v; i++) {
    for(const auto& [end_vert, weight] : g[i]) {
      if(std::find_if(g_undirected[end_vert].begin(), g_undirected[end_vert].end(),
                  [&](const auto& p) {return p.first == i;})
         == g_undirected[end_vert].end()) {
        g_undirected[end_vert].emplace_back(std::make_pair(i, weight));
      }
    }
  }

  std::vector<int> colour(v, 0);
  for (int i=0 ; i< v; i++) {
    if(colour[i] == 0) {
      if(!colourGraphUtil(g_undirected, i, n_c, colour)) {
        return {false, std::vector<int>()};
      }
    }
  }
  return {true, colour};
}

} // namespace grahAL

// using Pair = std::pair<int, double>;
//
// int main() {
// 	std::vector<std::vector<Pair>> g_v(5, std::vector<Pair>(3));
// 	g_v[0] = {std::make_pair(1,5), std::make_pair(4,4)};
// 	g_v[1] = {std::make_pair(3,7)};
// 	g_v[2] = {std::make_pair(1,8)};
// 	g_v[3] = {std::make_pair(4,9)};
// 	g_v[4] = {std::make_pair(0,5)};
//
// 	graphAL::GraphAL g_al(g_v);
// 	auto [found, path] = g_al.DFS(1,4);
// 	if(found) {
// 		std::cout << "Path: ";
// 		for(const auto& vert : path) {
// 			std::cout << vert << ' ';
// 		}
// 		std::cout << '\n';
// 	} else {
// 		std::cout << "path not found" << '\n';
// 	}
//     // Test for Prim
//   auto [found_prim, mst] = g_al.Prim();
//   if(found_prim) {
//     std::cout << "MST:" << '\n';
//     for(const auto& e : mst) {
//       std::cout << e.u << " --- " << e.v << '\n';
//     }
//   } else {
//     std::cout << "MST not found" << '\n';
//   }
//   auto mst = g_al.KruskalsAlgorithm();
//   std::cout << "MST:" << '\n';
//   for(const auto& e : mst) {
//     std::cout << e.u << " --- " << e.v << '\n';
//   }
//
//   auto cycle_found = g_al.UnionFindDetectCycle();
//   if(cycle_found) {
//     std::cout << "Cycle found" << '\n';
//   }else {
//     std::cout << "Cycle not found" << '\n';
//   }
//
//   auto flow = g_al.fordFulkerson(0,4);
//   std::cout << "Flow from source to sink: " << flow << '\n';
//
//   auto [mother_vertex_found, mv] = g_al.FindMotherVertex();
//
//   const auto[path_found, path_cost, path] = g_al.Dijkstra(0, 4);
//   const auto[path_found, p_md] = g_al.Dijkstra(0);
//   int n_c = 3;
//   const auto& [colourable, colour] = g_al.colourGraph(n_c);
//   std::cout << "The graph is colourable: " << colourable << " with " << n_c << " colours." << '\n';
//   for(const auto& c : colour) {
//     std::cout << c << '\n';
//   }
//   return 0;
// }
