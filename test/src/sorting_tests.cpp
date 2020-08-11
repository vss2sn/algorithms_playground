#include "algorithms/sort.hpp"
#include "algorithms/sort_iter.hpp"
#include "test/sorting_tests.hpp"

TEST_F(SortingTestRandomInt, InsertionSort) {
  sort::InsertionSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, InsertionSortWithCb) {
  sort::InsertionSortWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, InsertionSortWithRiCb) {
  sort::InsertionSortWithRiCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, InsertionSortWithBinarySearch) {
  sort::InsertionSortWithBinarySearch(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, InsertionSortWithBinarySearchWithCb) {
  sort::InsertionSortWithBinarySearchWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, InsertionSortOptimized) {
  sort::InsertionSortOptimized(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, InsertionSortIter) {
  sort::InsertionSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, InsertionSortOptimizedWithCb) {
  sort::InsertionSortOptimizedWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, BubbleSort) {
  sort::BubbleSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, BubbleSortIter) {
  sort::BubbleSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, MergeSort) {
  sort::MergeSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, QuickSort) {
  sort::QuickSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, SelectionSort) {
  sort::SelectionSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, SelectionSortIter) {
  sort::SelectionSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, BucketSort) {
  sort::BucketSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, CountingSort) {
  sort::CountingSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, HeapSort) {
  sort::HeapSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, STLSort) {
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
