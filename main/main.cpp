#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "algorithms/binary_search_tree.hpp"
#include "algorithms/binary_tree.hpp"
#include "algorithms/graphAM.hpp"
#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"

int main() {
	std::vector<int> to_sort{10,8,9,3,5,6,1,4,2,7};

	sort::SelectionSort(to_sort.begin(), to_sort.end());

	for(const auto& ele : to_sort) {
		std::cout << ele << ' ';
	}
	std::cout << '\n';
	return 0;
}
