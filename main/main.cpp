#include <iostream>

#include "algorithms/binary_search_tree.hpp"
#include "algorithms/binary_tree.hpp"
#include "algorithms/graphAL.hpp"
#include "algorithms/graphAM.hpp"
#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"

int main() {
	// std::vector<std::vector<double>> gv = {
	//     { 0, 1, 1, 0, 0},
	//     { 0, 0, 0, 0, 0},
	//     { 0, 1, 0, 0, 0},
	//     { 0, 1, 1, 0, 1},
	//     { 0, 1, 1, 0, 0},
	//   };
	// graphAM::GraphAM g(gv);

	std::vector<std::vector<std::pair<int, double>>> gv = {
			{ {1,1}, {2, 1} },
			{},
			{ {1,1} },
			{ {1,1}, {2, 1}, {4, 1} },
			{ {1,1}, {2, 1} },
		};
	graphAL::GraphAL g(gv);

	auto [sorted, order] = g.TopologicalSort();

	if(sorted) {
		std::cout << "Topological Sorted Order: ";
		for (const auto& ele : order) {
			std::cout << ele << ' ';
		}
		std::cout << '\n';
	} else {
		std::cout << "The graph is undirected/has a cycle" << '\n';
	}

	return 0;
}
