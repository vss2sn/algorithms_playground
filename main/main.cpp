#include <future>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "algorithms/binary_search_tree.hpp"
#include "algorithms/binary_tree.hpp"
#include "algorithms/graphAL.hpp"
#include "algorithms/graphAM.hpp"
#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"

using namespace std::chrono;

// int main() {
//  std::vector<int> to_sort{2,1,3,6,4,6,2,7,5,7,8,0};
// 	sort::MergeSortInPlace(to_sort.begin(), to_sort.end());
// 	assert(std::is_sorted(to_sort.begin(), to_sort.end()));
// 	return 0;
// }

// int main() {
//
// 	const int size_to_sort{100000};
//   std::vector<int> to_sort(size_to_sort);
//   std::random_device rd;
//
//   std::mt19937 gen(rd());
//   std::uniform_real_distribution<> dist(0, size_to_sort);
// 	int duration = 0;
//
// 	for(int i = 0; i < 10; i++) {
// 		std::cout << "Iteration: " << i << std::endl;
// 		std::generate(to_sort.begin(), to_sort.end(), [&](){return dist(gen);});
//
//     auto start = high_resolution_clock::now();
// 		// sort::MergeSortInPlace(to_sort);
// 		sort::MergeSortInPlace(to_sort.begin(), to_sort.end(), std::less<>());
// 		auto stop = high_resolution_clock::now();
//
// 		assert(std::is_sorted(to_sort.begin(), to_sort.end(), std::less<>()));
// 		duration += duration_cast<microseconds>(stop - start).count();
// 	}
//
// 	std::cout << "Time taken: " << duration << std::endl;
//
// 	return 0;
// }

using Pair = std::pair<int, double>;

int main() {
	std::vector<std::vector<Pair>> g_v(5, std::vector<Pair>());
	g_v[0] = {std::make_pair(2,1), std::make_pair(3,1)};
	g_v[1] = {std::make_pair(0,1)};
	g_v[2] = {std::make_pair(1,1)};
	g_v[3] = {std::make_pair(4,1)};
	g_v[4] = {};

	graphAL::GraphAL g_al(g_v);
  const auto components = g_al.KosarajuAlgorithm();
  for(auto& row : components) {
    std::cout << "Connected component: ";
    for(auto & element : row) {
      std::cout << element;
    }
    std::cout << '\n';
  }

  return 0;
}
