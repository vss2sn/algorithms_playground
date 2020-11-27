#include <iostream>

#include "algorithms/binary_search_tree.hpp"
#include "algorithms/binary_tree.hpp"
#include "algorithms/graphAL.hpp"
#include "algorithms/graphAM.hpp"
#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"

int main() {
	std::vector<std::vector<double>> gv = {
	    { 0, 1, 0, 0},
	    { 1, 0, 0, 0},
	    { 0, 0, 0, 1},
	    { 0, 0, 1, 0}
	  };
	graphAM::GraphAM g(gv);

	std::vector<std::pair<int, int>> ans =  g.FindBridges();
	std::cout << "Bridges:" << '\n';
	for(const auto& bridge : ans) {
		std::cout << bridge.first << ' ' << bridge.second << '\n';
	}

	return 0;
}
