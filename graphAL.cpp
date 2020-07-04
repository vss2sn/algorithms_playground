#include <iostream>
#include <queue>
#include <stack>
#include <vector>
#include <utility>

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
	auto [found, path] = g_al.DFS(1,4);
	if(found) {
		std::cout << "Path: ";
		for(const auto& vert : path) {
			std::cout << vert << ' ';
		}
		std::cout << '\n';
	} else {
		std::cout << "path not found" << '\n';
	}
}
