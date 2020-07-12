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
      if(parent[p.first] == -1) {
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

} // namespace grahAL


using Pair = std::pair<int, double>;

int main() {
	std::vector<std::vector<Pair>> g_v(5, std::vector<Pair>(3));
	g_v[0] = {std::make_pair(1,5), std::make_pair(2,1), std::make_pair(4,4)};
	g_v[1] = {std::make_pair(3,7), std::make_pair(4,2)};
	g_v[2] = {std::make_pair(1,8)};
	g_v[3] = {std::make_pair(4,9)};
	g_v[4] = {std::make_pair(1,1), std::make_pair(0,5)};

	graphAL::GraphAL g_al(g_v);
	// auto [found, path] = g_al.DFS(1,4);
	// if(found) {
	// 	std::cout << "Path: ";
	// 	for(const auto& vert : path) {
	// 		std::cout << vert << ' ';
	// 	}
	// 	std::cout << '\n';
	// } else {
	// 	std::cout << "path not found" << '\n';
	// }
  //   // Test for Prim
  // auto [found_prim, mst] = g_al.Prim();
  // if(found_prim) {
  //   std::cout << "MST:" << '\n';
  //   for(const auto& e : mst) {
  //     std::cout << e.u << " --- " << e.v << '\n';
  //   }
  // } else {
  //   std::cout << "MST not found" << '\n';
  // }
  auto mst = g_al.KruskalsAlgorithm();
  std::cout << "MST:" << '\n';
  for(const auto& e : mst) {
    std::cout << e.u << " --- " << e.v << '\n';
  }
}
