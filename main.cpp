#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "binary_tree.hpp"
#include "binary_search_tree.hpp"
#include "graphAM.hpp"
#include "sort.hpp"

int main() {
	std::vector<int> to_sort{10,8,9,3,5,6,1,4,2,7};
	// BubbleSort(to_sort);
	// InsertionSort(to_sort);
	// InsertionSortWithBinarySearch(to_sort);
	QuickSort(to_sort);
	for(const auto& ele : to_sort) {
		std::cout << ele << ' ';
	}
	std::cout << '\n';
	return 0;
}
