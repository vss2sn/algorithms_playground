#ifndef SORTING_TESTS_HPP
#define SORTING_TESTS_HPP

#include <algorithm>
#include <random>

#include <gtest/gtest.h>

class SortingTestRandomInt : public ::testing::Test {
protected:
  const int size_to_sort{10000};
  std::vector<int> to_sort_stl, to_sort_algo;
  std::random_device rd;

  SortingTestRandomInt()
  :to_sort_stl(size_to_sort,0){}

  void SetUp() override {
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, size_to_sort);
    std::generate(to_sort_stl.begin(), to_sort_stl.end(), [&](){return dist(gen);});
    to_sort_algo = to_sort_stl;
  }
};

#endif  // SORTING_TESTS_HPP
