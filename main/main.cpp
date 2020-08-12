#include <future>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "algorithms/binary_search_tree.hpp"
#include "algorithms/binary_tree.hpp"
#include "algorithms/graphAM.hpp"
#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"

using namespace std::chrono;

// int main() {
//   std::vector<int> to_sort{2,1,3,6,4,6,2,7,5,7,8,0};
// 	for(const auto p: to_sort) {
// 		std::cout << p << " ";
// 	}
// 	std::cout << '\n';
//
// 	sort::HeapSort(to_sort.begin(), to_sort.end());
// 	for(const auto p: to_sort) {
// 		std::cout << p << " ";
// 	}
// 	std::cout << '\n';
// 	assert(std::is_sorted(to_sort.begin(), to_sort.end()));
// 	return 0;
// }

int main() {

	const int size_to_sort{10000};
  std::vector<int> to_sort(size_to_sort);
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0, size_to_sort);
	int duration = 0;

	for(int i = 0; i < 1000; i++) {
		std::cout << "Iteration: " << i << std::endl;
		std::generate(to_sort.begin(), to_sort.end(), [&](){return dist(gen);});

		auto start = high_resolution_clock::now();
		// sort::QuickSort(to_sort);
		sort::QuickSort(to_sort.begin(), to_sort.end(), std::greater<>());
		auto stop = high_resolution_clock::now();

		assert(std::is_sorted(to_sort.begin(), to_sort.end(), std::greater<>()));
		duration += duration_cast<microseconds>(stop - start).count();
	}

	std::cout << "Time taken: " << duration << std::endl;

	return 0;
}