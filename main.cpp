#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "binary_tree.hpp"
#include "graphAM.hpp"

int main () {

  // Test for allPathsBetween
  graphAM::GraphAM g = {
    {0, 1, 1, 0, 0, 0},
    {1, 0, 1, 0, 0, 1},
    {1, 1, 0, 1, 1, 1},
    {0, 0, 1, 0, 1, 1},
    {0, 0, 1, 1, 0, 1},
    {0, 1, 1, 1, 1, 0}
  };

  auto [paths_found, paths] = graphAM::colourGraph(g, 4);
  return 0;
}
