#ifndef SORTING_TESTS_HPP
#define SORTING_TESTS_HPP

#include <algorithm>
#include <random>

#include <gtest/gtest.h>

class SortingTestRandomInt : public ::testing::Test {
protected:
  const int size_to_sort{100};
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

class SortingTestParamInt : public ::testing::TestWithParam<std::vector<int>>
{
public:
  std::vector<int> to_sort;

  void SetUp() override {
    to_sort = GetParam();
  }
};

class SortingTestParamChar : public ::testing::TestWithParam<std::vector<char>>
{
public:
  std::vector<char> to_sort;

  void SetUp() override {
    to_sort = GetParam();
  }
};

#endif  // SORTING_TESTS_HPP
