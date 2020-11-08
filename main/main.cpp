#include <iostream>

#include "algorithms/binary_search_tree.hpp"
#include "algorithms/binary_tree.hpp"
#include "algorithms/graphAL.hpp"
#include "algorithms/graphAM.hpp"
#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"

int main() {
	// std::vector<std::vector<double>> gv = {
	//     { 0, 1, 0, 0, 0, 0, 1},
	//     { 0, 0, 1, 0, 0, 0, 0},
	//     { 1, 0, 0, 1, 0, 0, 0},
	//     { 0, 0, 0, 0, 1, 0, 0},
	//     { 0, 0, 1, 0, 0, 1, 0},
	// 		{ 1, 0, 0, 0, 0, 0, 0},
	// 		{ 0, 0, 0, 1, 1, 0, 0},
	//   };
	// graphAM::GraphAM g(gv);

	std::vector<std::vector<std::pair<int, double>>> adj2(7);
  adj2[0].push_back(std::make_pair<int, double>(1, 1));
  adj2[0].push_back(std::make_pair<int, double>(6, 1));
  adj2[1].push_back(std::make_pair<int, double>(2, 1));
  adj2[2].push_back(std::make_pair<int, double>(0, 1));
  adj2[2].push_back(std::make_pair<int, double>(3, 1));
  adj2[3].push_back(std::make_pair<int, double>(4, 1));
  adj2[4].push_back(std::make_pair<int, double>(2, 1));
  adj2[4].push_back(std::make_pair<int, double>(5, 1));
  adj2[5].push_back(std::make_pair<int, double>(0, 1));
  adj2[6].push_back(std::make_pair<int, double>(4, 1));
  graphAL::GraphAL g(adj2);
	const auto path = g.HierholzersAlgorithm();
	for(const auto& ele : path) {
		std::cout << ele << ' ';
	}
	std::cout << '\n';
	return 0;
}
