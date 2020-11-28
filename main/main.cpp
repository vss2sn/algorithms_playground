#include <iostream>

#include "algorithms/binary_search_tree.hpp"
#include "algorithms/binary_tree.hpp"
#include "algorithms/graphAL.hpp"
#include "algorithms/graphAM.hpp"
#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"

int main() {
	// std::vector<std::vector<double>> gv = {
	//     { 0, 1, 0, 0},
	//     { 1, 0, 0, 1},
	//     { 0, 0, 0, 1},
	//     { 0, 1, 1, 0}
	//   };

	std::vector<std::vector<std::pair<int, double>>> gv = {
		{ {1,1}, {2,1}, {3,1} },
    { {0,1}, {2,1} },
    { {0,1}, {1,1} },
    { {0,1} },
  };
	graphAL::GraphAL g(gv);

	std::cout << g.IsBipartite() << '\n';

	return 0;
}
