#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "algorithms/graph.hpp"
#include "algorithms/graphAM.hpp"

TEST(GraphAM, BFS) {
  std::vector<std::vector<double>> g;
  std::vector<int> expected_result_v;
  std::tuple<bool, std::vector<int>> result;

  g = {
    { 1, 1, 1, 0, 0 },
    { 0, 1, 1, 0, 0 },
    { 0, 0, 1, 1, 0 },
    { 0, 0, 0, 1, 1 },
    { 0, 0, 0, 0, 1 },
  };

  auto g_am = graphAM::GraphAM(g);

  result = g_am.BFS(0, 4);
  ASSERT_EQ(std::get<0>(result), true);
  expected_result_v = {0, 2, 3, 4};
  ASSERT_EQ(std::get<1>(result), expected_result_v);

  result = g_am.BFS(4, 0);
  ASSERT_EQ(std::get<0>(result), false);

  result = g_am.BFS(2, 0);
  ASSERT_EQ(std::get<0>(result), false);

  result = g_am.BFS(0, 2);
  ASSERT_EQ(std::get<0>(result), true);
  expected_result_v = {0, 2};
  ASSERT_EQ(std::get<1>(result), expected_result_v);

  result = g_am.BFS(2, 4);
  ASSERT_EQ(std::get<0>(result), true);
  expected_result_v = {2, 3, 4};
  ASSERT_EQ(std::get<1>(result), expected_result_v);
}

TEST(GraphAM, DFS) {
  std::vector<std::vector<double>> g;
  std::vector<int> expected_result_v;
  std::tuple<bool, std::vector<int>> result;

  g = {
    { 1, 1, 1, 0, 0 },
    { 0, 1, 1, 0, 0 },
    { 0, 0, 1, 1, 0 },
    { 0, 0, 0, 1, 1 },
    { 0, 0, 0, 0, 1 },
  };

  auto g_am = graphAM::GraphAM(g);

  result = g_am.DFS(0, 4);
  expected_result_v = {0, 2, 3, 4};

  ASSERT_EQ(std::get<0>(result), true);
  ASSERT_EQ(std::get<1>(result), expected_result_v);

  result = g_am.DFS(4, 0);
  ASSERT_EQ(std::get<0>(result), false);

  result = g_am.DFS(2, 0);
  ASSERT_EQ(std::get<0>(result), false);

  result = g_am.DFS(0, 2);
  ASSERT_EQ(std::get<0>(result), true);
  expected_result_v = {0, 2};
  ASSERT_EQ(std::get<1>(result), expected_result_v);

  result = g_am.DFS(2, 4);
  ASSERT_EQ(std::get<0>(result), true);
  expected_result_v = {2, 3, 4};
  ASSERT_EQ(std::get<1>(result), expected_result_v);
}

TEST(GraphAM, Prim1) {
  std::vector<std::vector<double>> g;
  std::tuple<bool, std::vector<Edge>> expected_result;
  std::tuple<bool, std::vector<Edge>> result;

  g = {
    { 1, 1, 1, 0, 0 },
    { 0, 1, 1, 0, 0 },
    { 0, 0, 1, 1, 0 },
    { 0, 0, 0, 1, 1 },
    { 0, 0, 0, 0, 1 },
  };
  auto g_am = graphAM::GraphAM(g);

  result = g_am.Prim();
  expected_result = {false, {}};
  ASSERT_EQ(std::get<0>(result), std::get<0>(expected_result));
}

TEST(GraphAM, Prim2) {
  std::vector<std::vector<double>> g;
  std::tuple<bool, std::vector<Edge>> expected_result;
  std::tuple<bool, std::vector<Edge>> result;

  g = {
    { 1, 1, 1, 0, 0 },
    { 1, 1, 1, 0, 0 },
    { 1, 1, 1, 1, 0 },
    { 0, 0, 1, 1, 1 },
    { 0, 0, 0, 1, 1 },
  };
  auto g_am = graphAM::GraphAM(g);

  result = g_am.Prim();
  expected_result = {true, {{0,1,1},{0,2,1},{2,3,1},{3,4,1}}};
  ASSERT_EQ(std::get<0>(result), std::get<0>(expected_result));
  ASSERT_EQ(std::get<1>(result), std::get<1>(expected_result));
}

TEST(GraphAM, FordFulkerson) {
  std::vector<std::vector<double>> g;
  std::vector<int> expected_result_v;
  std::tuple<bool, std::vector<int>> result;

  g = {
    {0, 16, 13, 0, 0, 0},
    {0, 0, 10, 12, 0, 0},
    {0, 4, 0, 0, 14, 0},
    {0, 0, 9, 0, 0, 20},
    {0, 0, 0, 7, 0, 4},
    {0, 0, 0, 0, 0, 0}
  };

  auto g_am = graphAM::GraphAM(g);
	double flow = g_am.fordFulkerson(0, 5);
  ASSERT_EQ(flow, 23);
}

TEST(GraphAM, MotherVertex) {
  std::vector<std::vector<double>> g;
  std::tuple<bool, int> expected_result;
  std::tuple<bool, int> result;

  // Test for mother vertex
  g = {
    {0, 1, 1, 0, 0, 0, 0},
    {0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0, 0},
    {0, 0, 1, 0, 0, 0, 1},
    {1, 0, 0, 0, 1, 0, 0}
  };

  auto g_am = graphAM::GraphAM(g);
  expected_result = {true, 5};

  result = g_am.FindMotherVertex();
  ASSERT_EQ(std::get<0>(result), std::get<0>(expected_result));
  ASSERT_EQ(std::get<1>(result), std::get<1>(expected_result));
}

TEST(GraphAM, UnionFind) {
  std::vector<std::vector<double>> g;
  bool expected_result = false;
  bool result = false;

  // Test for UnionFind and UnionFindRC
  g = {
    {0, 1, 1, 0, 0, 0},
    {0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0}
  };

  auto g_am = graphAM::GraphAM(g);

  expected_result = true;
  result = g_am.UnionFindDetectCycle();
  ASSERT_EQ(result, expected_result);
  result = g_am.UnionFindRCDetectCycle();
  ASSERT_EQ(result, expected_result);
}

TEST(GraphAM, AllPathsBetween) {
  std::vector<std::vector<double>> g;
  std::tuple<bool, std::vector<std::vector<int>>> expected_result;
  std::tuple<bool, std::vector<std::vector<int>>> result;

  // Test for allPathsBetween
  g = {
    {0, 1, 1, 0, 0, 0},
    {1, 0, 1, 0, 0, 1},
    {1, 1, 0, 1, 1, 1},
    {0, 0, 1, 0, 1, 1},
    {0, 0, 1, 1, 0, 1},
    {0, 1, 1, 1, 1, 0}
  };

  auto g_am = graphAM::GraphAM(g);
  expected_result = {true, {{5, 1, 0},
                            {5, 1, 2, 0},
                            {5, 2, 0},
                            {5, 2, 1, 0},
                            {5, 3, 2, 0},
                            {5, 3, 2, 1, 0},
                            {5, 3, 4, 2, 0},
                            {5, 3, 4, 2, 1, 0},
                            {5, 4, 2, 0},
                            {5, 4, 2, 1, 0},
                            {5, 4, 3, 2, 0},
                            {5, 4, 3, 2, 1, 0 }}};
  result = g_am.allPathsBetween(0, 5);
  ASSERT_EQ(std::get<0>(result), std::get<0>(expected_result));
  ASSERT_EQ(std::get<1>(result), std::get<1>(expected_result));
}

TEST(GraphAM, ColourGraph) {
  std::vector<std::vector<double>> g;
  std::tuple<bool, std::vector<int>> expected_result;
  std::tuple<bool, std::vector<int>> result;

  // Test for colourGraph
  g = {
    {0, 1, 1, 0, 0, 0},
    {1, 0, 1, 0, 0, 1},
    {1, 1, 0, 1, 1, 1},
    {0, 0, 1, 0, 1, 1},
    {0, 0, 1, 1, 0, 1},
    {0, 1, 1, 1, 1, 0}
  };

  auto g_am = graphAM::GraphAM(g);
  expected_result = {true, {1, 2, 3, 1, 2, 4}};
  result = g_am.colourGraph(4);
  ASSERT_EQ(std::get<0>(result), std::get<0>(expected_result));
  ASSERT_EQ(std::get<1>(result), std::get<1>(expected_result));
}

TEST(GraphAM, Dijkstra) {
  std::vector<std::vector<double>> g;
  std::tuple<bool, double, std::vector<int>> expected_result;
  std::tuple<bool, double, std::vector<int>> result;

g = {
   {0, 1, 1, 0, 0, 0},
   {1, 0, 1, 0, 0, 0},
   {1, 1, 0, 1, 1, 0},
   {0, 0, 1, 0, 1, 0},
   {0, 0, 1, 1, 0, 1},
   {0, 1, 1, 1, 1, 0}
  };
  auto g_am = graphAM::GraphAM(g);
  expected_result = {true, 3, {0, 2, 4, 5}};
  result = g_am.Dijkstra(0, 5);

  ASSERT_EQ(std::get<0>(result), std::get<0>(expected_result));
  ASSERT_EQ(std::get<1>(result), std::get<1>(expected_result));
  ASSERT_EQ(std::get<2>(result), std::get<2>(expected_result));
}

TEST(GraphAM, DijkstraAll) {
  std::vector<std::vector<double>> g;
  std::tuple<bool, std::vector<std::pair<int, double>>> expected_result;
  std::tuple<bool, std::vector<std::pair<int, double>>> result;

  g = {
     {0, 1, 1, 0, 0, 0},
     {1, 0, 1, 0, 0, 0},
     {1, 1, 0, 1, 1, 0},
     {0, 0, 1, 0, 1, 0},
     {0, 0, 1, 1, 0, 0},
     {0, 1, 1, 1, 0, 0}
    };
  auto g_am = graphAM::GraphAM(g);
  expected_result = {true, {{1, 2}, {0, 1}, {0, 1}, {2, 2}, {2, 2},  {-1, std::numeric_limits<double>::max()} }};
  result = g_am.Dijkstra(0);
  ASSERT_EQ(std::get<0>(result), std::get<0>(expected_result));
  ASSERT_EQ(std::get<1>(result), std::get<1>(expected_result));
}

TEST(GraphAM, Kosaraju) {
  std::vector<std::vector<double>> g;
  std::vector<std::vector<int>> expected_result;
  std::vector<std::vector<int>> result;

  g = {
   {0, 0, 1, 1, 0},
   {1, 0, 0, 0, 0},
   {0, 1, 0, 0, 0},
   {0, 0, 0, 0, 1},
   {0, 0, 0, 0, 0},
  };
  auto g_am = graphAM::GraphAM(g);
  result = g_am.KosarajuAlgorithm();
  expected_result = {{0, 1, 2}, {3}, {4}};
  ASSERT_EQ(result, expected_result);
}

TEST(GraphAM, StronglyConnectedKosaraju) {
  std::vector<std::vector<double>> g;
  bool expected_result;
  bool result;

  g = {
   {0, 0, 1, 1, 0},
   {1, 0, 0, 0, 0},
   {1, 1, 0, 0, 0},
   {1, 0, 0, 0, 1},
   {0, 0, 0, 1, 0},
  };
  auto g_am = graphAM::GraphAM(g);
  expected_result = true;
  result = g_am.StronglyConnectedKosaraju();
  ASSERT_EQ(result, expected_result);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
