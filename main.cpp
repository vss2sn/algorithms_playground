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
	// sort::BubbleSort(to_sort);

	// sort::InsertionSort(to_sort);
  // sort::InsertionSortWithCb(to_sort);
  // sort::InsertionSortWithRiCb(to_sort);
	// sort::InsertionSortWithBinarySearch(to_sort);
  // sort::InsertionSortWithBinarySearchWithCb(to_sort);

  // sort::MergeSort(to_sort);

  // sort::QuickSort(to_sort);

  // sort::SelectionSort(to_sort);

	// sort::BucketSort(to_sort);

	sort::CountingSort(to_sort);

	for(const auto& ele : to_sort) {
		std::cout << ele << ' ';
	}
	std::cout << '\n';
	return 0;
}
