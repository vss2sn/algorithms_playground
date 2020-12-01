#include <iostream>
#include <queue>
#include <stack>
#include <utility>
#include <unordered_set>
#include <vector>

#include <thread>
#include <chrono>

#include "algorithms/graphAL.hpp"

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
    pq.push({source, end_vert, weight});
  }

  std::vector<Edge> mst;
  while(!pq.empty() && mst.size() < v) {
    Edge current_edge = pq.top();
    pq.pop();
    if(parent[current_edge.v] == -1 ) {
      parent[current_edge.v] = current_edge.u;
      mst.push_back(current_edge);
      for(const auto& [end_vert, weight] : g[current_edge.v]) {
        pq.push({current_edge.v, end_vert, weight});
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
      edge_pq.push({start_vert, end_vert, weight});
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

bool GraphAL::BellmanFord(const int source) const {
  std::vector<double> distances(v, std::numeric_limits<double>::max());
  std::vector<double> predecessors(v, -1);

  distances[source] = 0;

  for(int k = 0; k < v; ++k) {
    for(int i = 0; i < v; ++i) {
      for(const auto& [end_vert, weight] : g[i]) {
        if(double delta = distances[i] + weight; distances[end_vert] > delta) {
          distances[end_vert] = delta;
          predecessors[end_vert] = i;
        }
      }
    }
  }

  for(int i = 0; i < v; ++i) {
    for(const auto& [end_vert, weight] : g[i]) {
      if(distances[end_vert] > distances[i] + weight) {
        std::cout << __FUNCTION__ << __LINE__ << "Negative cycle detected" << '\n';
        return true;
      }
    }
  }

  return true;
}

std::vector<std::vector<std::pair<int, double>>> GraphAL::invertGraph() const {
  std::vector<std::vector<std::pair<int, double>>>  g_inv(v, std::vector<std::pair<int, double>>());
  for(int i = 0; i < v; ++i) {
    for(const auto& [end_vert, weight] : g[i]) {
      g_inv[end_vert].push_back(std::make_pair(i, weight));
    }
  }
  return g_inv;
}

void GraphAL::DFSUtilWithFinishTime(int source, std::vector<bool>& visited, std::stack<int>& visit_order) const{
  visited[source] = true;
  for(const auto& [end_vert, weights] : g[source]) {
    if(!visited[end_vert]) {
      DFSUtilWithFinishTime(end_vert, visited, visit_order);
    }
  }
  visit_order.push(source);
}

std::vector<std::vector<int>> GraphAL::KosarajuAlgorithmUtil() {
  std::vector<bool> visited(v, false);
  std::stack<int> visit_order;
  for(int i = 0; i < v; ++i) {
    if(!visited[i]) {
      DFSUtilWithFinishTime(i, visited, visit_order);
    }
  }

  g = invertGraph();
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

std::vector<std::vector<int>> GraphAL::KosarajuAlgorithm() const {
  GraphAL g_al(g);
  return g_al.KosarajuAlgorithmUtil();
}

bool GraphAL::StronglyConnectedKosarajuUtil() {
  std::vector<bool> visited(v, false);
  std::stack<int> visit_order;
  DFSUtil(0, visited);

  for(const auto& ele : visited) {
    if(!ele) return false;
  }

  g = invertGraph();
  std::fill(visited.begin(), visited.end(), false);
  DFSUtil(0, visited);

  for(const auto& ele : visited) {
    if(!ele) return false;
  }

  return true;
}

bool GraphAL::StronglyConnectedKosaraju() const {
  GraphAL g_al(g);
  return g_al.StronglyConnectedKosarajuUtil();
}


std::vector<int> GraphAL::HierholzersAlgorithmUtil() {
  // check degrees
  std::vector<int> indeg(v, 0), outdeg(v, 0);
  for(int i = 0 ; i < v; i++) {
    for(const auto& edge : g[i] ) {
        ++indeg[edge.first];
        ++outdeg[i];
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
      g[v1].emplace_back(v2, 1);
    } else {
      g[v2].emplace_back(v1, 1);
    }
  }

  std::stack<int> vertices; // vertices to expand
  std::vector<int> eulerian_path;
  vertices.push(first);
  while(!vertices.empty()) {
    const int edge_start = vertices.top();
    auto iter = std::find_if_not(std::begin(g[edge_start]), std::end(g[edge_start]),
      [](auto p) { return p.second == 0; });
    if(iter == g[edge_start].end()) {
      eulerian_path.push_back(edge_start);
      vertices.pop();
    } else {
      iter->second = 0;
      vertices.push(iter->first);
    }
  }

  // remove added edge if added
    const int s1 = eulerian_path.size() - 1;
    if(v1 != -1) {
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
      if(ele.second != 0) {
        std::cout << "There is an edge left" << '\n';
        return std::vector<int>();
      }
    }
  }

  std::reverse(eulerian_path.begin(), eulerian_path.end());
  return eulerian_path;
}

std::vector<int> GraphAL::HierholzersAlgorithm() const {
  GraphAL mod(g);
  return mod.HierholzersAlgorithmUtil();
}

void GraphAL::ArticulationPointsUtil(const int vert1, std::vector<bool>& visited, std::vector<int>& parent,
  std::vector<int>& tod, std::vector<int>& low, std::unordered_set<int>& articulation_points, int time) const {
  int children = 0;
  visited[vert1] = true;
  low[vert1] = tod[vert1] = ++time;
  const int n_edges = g[vert1].size();
  for(int j = 0; j < n_edges; j++) {
    const auto [vert2, edge_weight] = g[vert1][j];
    if(!visited[j]) {
      ++children;
      parent[vert2] = vert1;
      ArticulationPointsUtil(vert2, visited, parent, tod, low, articulation_points, time);
      low[vert1] = std::min(low[vert1], low[vert2]);
      if(low[vert1] >= tod[vert2] && parent[vert1] != -1) {
        articulation_points.insert(vert1);
      } else {
        low[vert1] = std::min(low[vert1], tod[vert2]);
      }
    }
    if(parent[vert1] != -1 && children > 1) {
      articulation_points.insert(vert1);
    }
  }
}

std::unordered_set<int> GraphAL::ArticulationPoints() const {

  // While this algorithm will correctly find all teh articulation points as long as
  // for each edge between 2 vertices there exists an edge runnin in the opposite direction
  // regardless of weight, this algorithm returns when the graph is not symetric to emphasize
  // the fact that this algorithm is meant to be run n undirected graphs
  for (int vert1 = 0; vert1 < v; ++vert1) {
    for(int j = 0; j <  g[vert1].size(); ++j) {
      const auto [vert2, edge_weight2] = g[vert1][j];
      bool found = false;
      for(const auto& [vert3, edge_weight3] : g[vert2]) {
        if(vert3 == vert1 && edge_weight3 == edge_weight2) {
          found= true;
          break;
        }
      }
      if(found == false) {
        std::cout << "Edge from " << vert1 << " to " << vert2 <<
        " not symetric. The algorithm will fail. " << '\n';
        return std::unordered_set<int>();
      }
    }
  }

  std::vector<bool> visited(v, false);
  std::vector<int> parent(v, -1);
  std::vector<int> tod(v, std::numeric_limits<int>::max());
  std::vector<int> low(v, std::numeric_limits<int>::max());
  std::unordered_set<int> articulation_points;
  int time = 0;

  for(int i = 0; i < v; i++) {
    if(!visited[i]) {
      ArticulationPointsUtil(i, visited, parent, tod, low, articulation_points, time);
    }
  }

  return articulation_points;
}

bool GraphAL::TopologicalSortUtil(std::vector<int>& sorted,
  std::vector<bool>& visited,
  std::vector<bool>& this_cycle,
  const int vert) const {

  this_cycle[vert] = true;
  visited[vert] = true;
  for(const auto& [vert2, weight] : g[vert]) {
    if(this_cycle[vert2] && visited[vert2]) {
      return false;
    } else if(!visited[vert2]) {
      if(!TopologicalSortUtil(sorted, visited, this_cycle, vert2)) {
        return false;
      }
    }
  }
  this_cycle[vert] = false;
  sorted.push_back(vert);
  return true;
}

std::tuple<bool, std::vector<int>> GraphAL::TopologicalSort() const {
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

bool GraphAL::IsBipartite() const {
  std::vector<int> side(v, -1);
  std::queue<std::pair<int, double>> q;

  for(int vert1 = 0; vert1 < v; vert1++) {
    if(side[vert1] == -1) {
      q.push({vert1, 0});
      side[vert1] = 0;
      while(!q.empty()) {
        const int vert2 = q.front().first;
        q.pop();
        for (const auto& [vert3, weight] : g[vert2]) {
          if(side[vert3] == -1) {
            side[vert3] = side[vert2] ^ 1; // xor
            q.push({vert3, 0});
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

bool GraphAL::DivideIntoTwoCliques() const {
  std::vector<std::vector<std::pair<int, double>>> g2(v);
  for(int row = 0; row < v; ++row) {
    int cur = 0;
    for(int edge = 0; edge < g[row].size(); ++edge) {
      if(g[row][edge].first == row) {
        std::cout << "Self edge exists" << '\n';
        return false;
      }
      while(cur < g[row][edge].first) {
        if(cur != row) {
          g2[row].push_back({cur, 1});
        }
        ++cur;
      }
    }
  }

  GraphAL gcomp(g2);
  return gcomp.IsBipartite();
}

bool GraphAL::HamiltonianPathBacktrackUtil(const int vert, std::vector<bool>& visited, int level, std::vector<int>& path) const {
  ++level;
  visited[vert] = true;
  path.push_back(vert);
  if(level == v) return true;
  for(const auto& [vert2, weight] : g[vert]) {
    if(!visited[vert2] && HamiltonianPathBacktrackUtil(vert2, visited, level, path) && vert2 != vert) {
      return true;
    }
  }
  --level;
  visited[vert] = false;
  path.pop_back();
  return false;
}

std::tuple<bool, std::vector<int>> GraphAL::HamiltonianPath() const {
  std::vector<bool> visited(v, false);
  std::vector<int> path;
  path.reserve(v);
  int level = 0;
  bool path_found = false;
  for(int vert = 0; vert < v; ++vert) {
    path_found = HamiltonianPathBacktrackUtil(vert, visited, level, path);
    if(path_found) {
      break;
    }
  }
  if(path_found) {
    std::cout << "Path: ";
    for(const auto& ele : path) {
      std::cout << ele << ' ';
    }
    std::cout << '\n';
  } else {
    std::cout << "No hamiltonian path found" << '\n';
  }
  return {path_found, path};
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
