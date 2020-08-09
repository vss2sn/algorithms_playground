#ifndef SORTING_TESTS_HPP
#define SORTING_TESTS_HPP

#include <algorithm>
#include <random>

#include <gtest/gtest.h>

class SortingTestRandomInt : public ::testing::Test {
protected:
  const int size_to_sort{10000};
  std::vector<int> to_sort;
  std::random_device rd;

  SortingTestRandomInt()
  :to_sort(size_to_sort,0){}

  void SetUp() override {
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, size_to_sort);
    std::generate(to_sort.begin(), to_sort.end(), [&](){return dist(gen);});
  }
};

#endif  // SORTING_TESTS_HPP
