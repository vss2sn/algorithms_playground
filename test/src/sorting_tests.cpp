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

TEST_F(SortingTestRandomInt, MergeSortIter) {
  sort::MergeSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, MergeSortInPlace) {
  sort::MergeSortInPlace(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, MergeSortInPlaceIter) {
  sort::MergeSortInPlace(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, QuickSort) {
  sort::QuickSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, QuickSortIter) {
  sort::QuickSort(to_sort.begin(), to_sort.end());
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

TEST_F(SortingTestRandomInt, HeapSortIter) {
  sort::HeapSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_F(SortingTestRandomInt, STLSort) {
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSort) {
  sort::InsertionSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSortWithCb) {
  sort::InsertionSortWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSortWithRiCb) {
  sort::InsertionSortWithRiCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSortWithBinarySearch) {
  sort::InsertionSortWithBinarySearch(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSortWithBinarySearchWithCb) {
  sort::InsertionSortWithBinarySearchWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSortOptimized) {
  sort::InsertionSortOptimized(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSortIter) {
  sort::InsertionSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, InsertionSortOptimizedWithCb) {
  sort::InsertionSortOptimizedWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, BubbleSort) {
  sort::BubbleSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, BubbleSortIter) {
  sort::BubbleSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, MergeSort) {
  sort::MergeSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, MergeSortIter) {
  sort::MergeSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, MergeSortInPlace) {
  sort::MergeSortInPlace(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, MergeSortInPlaceIter) {
  sort::MergeSortInPlace(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, QuickSort) {
  sort::QuickSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, QuickSortIter) {
  sort::QuickSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, SelectionSort) {
  sort::SelectionSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, SelectionSortIter) {
  sort::SelectionSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, BucketSort) {
  sort::BucketSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, CountingSort) {
  sort::CountingSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, HeapSort) {
  sort::HeapSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, HeapSortIter) {
  sort::HeapSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamInt, STLSort) {
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

INSTANTIATE_TEST_CASE_P(SortTest, SortingTestParamInt,
  ::testing::Values(std::vector<int>{0, 0, 0, 0, 0},
                    std::vector<int>{5, 4, 3, 2, 1},
                    std::vector<int>{1, 2, 3, 4, 5},
                    std::vector<int>{1, 2, 3, 2, 1},
                    std::vector<int>{1, 2, 3, 5, 4},
                    std::vector<int>{1, 2, 3, 4, 5, 6},
                    std::vector<int>{6, 5, 4, 3, 2, 1},
                    std::vector<int>{1, 2, 3, 3, 2, 1},
                    std::vector<int>{1, 2, 3, 4, 6, 5}));

TEST_P(SortingTestParamChar, InsertionSort) {
  sort::InsertionSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, InsertionSortWithCb) {
  sort::InsertionSortWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, InsertionSortWithRiCb) {
  sort::InsertionSortWithRiCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, InsertionSortWithBinarySearch) {
  sort::InsertionSortWithBinarySearch(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, InsertionSortWithBinarySearchWithCb) {
  sort::InsertionSortWithBinarySearchWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, InsertionSortOptimized) {
  sort::InsertionSortOptimized(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, InsertionSortIter) {
  sort::InsertionSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, InsertionSortOptimizedWithCb) {
  sort::InsertionSortOptimizedWithCb(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, BubbleSort) {
  sort::BubbleSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, BubbleSortIter) {
  sort::BubbleSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, MergeSort) {
  sort::MergeSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, MergeSortIter) {
  sort::MergeSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, MergeSortInPlace) {
  sort::MergeSortInPlace(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, MergeSortInPlaceIter) {
  sort::MergeSortInPlace(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, QuickSort) {
  sort::QuickSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, QuickSortIter) {
  sort::QuickSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, SelectionSort) {
  sort::SelectionSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, SelectionSortIter) {
  sort::SelectionSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, BucketSort) {
  sort::BucketSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, CountingSort) {
  sort::CountingSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, HeapSort) {
  sort::HeapSort(to_sort);
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, HeapSortIter) {
  sort::HeapSort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

TEST_P(SortingTestParamChar, STLSort) {
  std::sort(to_sort.begin(), to_sort.end());
  ASSERT_TRUE(std::is_sorted(to_sort.begin(), to_sort.end()));
}

INSTANTIATE_TEST_CASE_P(SortTest, SortingTestParamChar,
  ::testing::Values(std::vector<char>{'a', 'b', 'c', 'd', 'e'},
                    std::vector<char>{'e', 'd', 'c', 'b', 'a'},
                    std::vector<char>{'a', 'a', 'a', 'a', 'a'},
                    std::vector<char>{'a', 'b', 'c', 'b', 'a'}
                  ));

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
