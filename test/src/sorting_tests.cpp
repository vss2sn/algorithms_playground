#include "algorithms/sort.hpp"
#include "test/sorting_tests.hpp"

TEST_F(SortingTestRandomInt, InsertionSort) {
  sort::InsertionSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, InsertionSortWithCb) {
  sort::InsertionSortWithCb(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, InsertionSortWithRiCb) {
  sort::InsertionSortWithRiCb(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, InsertionSortWithBinarySearch) {
  sort::InsertionSortWithBinarySearch(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, InsertionSortWithBinarySearchWithCb) {
  sort::InsertionSortWithBinarySearchWithCb(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, InsertionSortOptimized) {
  sort::InsertionSortOptimized(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, InsertionSortOptimizedWithCb) {
  sort::InsertionSortOptimizedWithCb(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, BubbleSort) {
  sort::BubbleSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, MergeSortInPlace) {
  sort::MergeSortInPlace(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, MergeSort) {
  sort::MergeSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, QuickSort) {
  sort::QuickSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, SelectionSort) {
  sort::SelectionSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, BucketSort) {
  sort::BucketSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, CountingSort) {
  sort::CountingSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

TEST_F(SortingTestRandomInt, HeapSort) {
  sort::HeapSort(to_sort_algo);
  std::sort(to_sort_stl.begin(), to_sort_stl.end());
  ASSERT_EQ(to_sort_algo, to_sort_stl);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
