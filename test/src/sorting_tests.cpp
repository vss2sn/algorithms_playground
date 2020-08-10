#include "algorithms/sort.hpp"
#include "test/sorting_tests.hpp"

TEST_F(SortingTestRandomInt, InsertionSort) {
  auto to_sort_insertion_sort(to_sort);
  sort::InsertionSort(to_sort_insertion_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(to_sort_insertion_sort, to_sort);
}

TEST_F(SortingTestRandomInt, InsertionSortWithCb) {
  auto to_sort_insertion_sort_with_cb(to_sort);
  sort::InsertionSortWithCb(to_sort_insertion_sort_with_cb);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(to_sort_insertion_sort_with_cb, to_sort);
}

TEST_F(SortingTestRandomInt, InsertionSortWithRiCb) {
  auto to_sort_insertion_sort_with_ri_cb(to_sort);
  sort::InsertionSortWithRiCb(to_sort_insertion_sort_with_ri_cb);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(to_sort_insertion_sort_with_ri_cb, to_sort);
}

TEST_F(SortingTestRandomInt, InsertionSortWithBinarySearch) {
  auto to_sort_insertion_sort_with_binary_search(to_sort);
  sort::InsertionSortWithBinarySearch(to_sort_insertion_sort_with_binary_search);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(to_sort_insertion_sort_with_binary_search, to_sort);
}

TEST_F(SortingTestRandomInt, InsertionSortWithBinarySearchWithCb) {
  auto to_sort_insertion_sort_with_binary_search_with_cb(to_sort);
  sort::InsertionSortWithBinarySearchWithCb(to_sort_insertion_sort_with_binary_search_with_cb);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(to_sort_insertion_sort_with_binary_search_with_cb, to_sort);
}

TEST_F(SortingTestRandomInt, InsertionSortOptimized) {
  auto to_sort_insertion_sort_optimized(to_sort);
  sort::InsertionSortOptimized(to_sort_insertion_sort_optimized);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(to_sort_insertion_sort_optimized, to_sort);
}

TEST_F(SortingTestRandomInt, InsertionSortOptimizedWithCb) {
  auto to_sort_insertion_sort_optimized_with_cb(to_sort);
  sort::InsertionSortOptimizedWithCb(to_sort_insertion_sort_optimized_with_cb);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(to_sort_insertion_sort_optimized_with_cb, to_sort);
}

TEST_F(SortingTestRandomInt, BubbleSort) {
  auto bubble_sort(to_sort);
  sort::BubbleSort(bubble_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(bubble_sort, to_sort);
}

TEST_F(SortingTestRandomInt, MergeSort) {
  auto merge_sort(to_sort);
  sort::MergeSort(merge_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(merge_sort, to_sort);
}

TEST_F(SortingTestRandomInt, QuickSort) {
  auto quick_sort(to_sort);
  sort::QuickSort(quick_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(quick_sort, to_sort);
}

TEST_F(SortingTestRandomInt, SelectionSort) {
  auto selection_sort(to_sort);
  sort::SelectionSort(selection_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(selection_sort, to_sort);
}

TEST_F(SortingTestRandomInt, BucketSort) {
  auto bucket_sort(to_sort);
  sort::BucketSort(bucket_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(bucket_sort, to_sort);
}

TEST_F(SortingTestRandomInt, CountingSort) {
  auto counting_sort(to_sort);
  sort::CountingSort(counting_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(counting_sort, to_sort);
}

TEST_F(SortingTestRandomInt, HeapSort) {
  auto heap_sort(to_sort);
  sort::HeapSort(heap_sort);
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_EQ(heap_sort, to_sort);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
