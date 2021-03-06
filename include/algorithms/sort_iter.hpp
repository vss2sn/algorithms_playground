#include <algorithm>

namespace sort {

template<typename It, typename Compare = std::less<>>
void BubbleSort(It begin, It end, Compare compare = Compare{}) {
  for(auto it_out = begin; it_out != end; ++it_out) {
    for(auto it_in = begin; it_in < it_out; ++it_in) {
      if(compare(*it_out, *it_in)) {
        std::iter_swap(it_in, it_out);
      }
    }
  }
}

template<typename It, typename Compare = std::less<>>
void InsertionSort(It begin, It end, Compare compare = Compare{}) {
  for(auto it = begin; it != end; ++it) {
    std::rotate(std::upper_bound(begin, it, *it, compare), it, std::next(it));
  }
}

template<typename It, typename Compare = std::less<>>
void SelectionSort(It begin, It end, Compare compare = Compare{}) {
  for(auto it = begin; it != end; ++it) {
    std::iter_swap(it, std::min_element(it, end, compare));
  }
}

// template<typename It, typename Compare = std::less<>>
// void MergeSortMergeUtil(It begin, It mid, It end, Compare compare = Compare{}) {
//   std::vector<typename It::value_type> temp;
//   auto dist = (std::distance(begin, end));
//   temp.reserve(dist);
//   It it_1 = begin, it_2 = mid;
//   while(it_1 != mid && it_2 != end) {
//     if(compare(*it_1, *it_2)) {
//       temp.push_back(*it_1++);
//     } else {
//       temp.push_back(*it_2++);
//     }
//   }
//   temp.insert(temp.end(), it_1, mid);
//   temp.insert(temp.end(), it_2, end);
//   std::copy(temp.begin(), temp.end(), begin);
// }

template<typename It, typename Compare = std::less<>>
void MergeSortMergeUtil(It begin, It mid, It end, Compare compare = Compare{}) {

  std::vector<typename It::value_type> left(begin, mid);
	std::vector<typename It::value_type> right(mid, end);
  It current = begin;
  It it_1 = left.begin(), it_2 = right.begin();
  while(it_1 != left.end() && it_2 != right.end()) {
    if(compare(*it_1, *it_2)) {
      *current++ = *it_1++;
    } else {
      *current++ = *it_2++;
    }
  }
  auto temp = std::copy(it_1, left.end(), current);
  std::copy(it_2, right.end(), temp);
}

template<typename It, typename Compare = std::less<>>
void MergeSort(It begin, It end, Compare compare = Compare{}) {
  auto dist = std::distance(begin, end);
  if(dist < 2) {
    return;
  }
  auto mid = std::next(begin, dist/2);
  MergeSort(begin, mid, compare);
  MergeSort(mid, end, compare);
  MergeSortMergeUtil(begin, mid, end, compare);
}

// template<typename It, typename Compare = std::less<>>
// void MergeSortInPlaceUtil(It begin, It mid, It end, Compare compare = Compare{}){
//   // auto begin_1 = begin;
//   // auto begin_2 = mid;
//   while(begin != mid && mid != end) {
//     if(compare(*begin, *mid)) {
//       ++begin;
//     } else {
//       // std::rotate(begin_1++, begin_2++, begin_2);
//       std::copy_backward(begin++, mid++, mid);
//       // std::advance(begin_1, 1);
//       // std::advance(begin_2, 1);
//     }
//   }
// }


template<typename It, typename Compare = std::less<>>
void MergeSortInPlaceUtil(It begin, It mid, It end, Compare compare = Compare{}){
  while(begin != mid && mid != end) {
    if(compare(*begin, *mid)) {
      ++begin;
    } else {
      // std::rotate(begin++, mid, std::next(mid, 1));
      std::copy_backward(begin++, mid, std::next(mid, 1));
      ++mid;
    }
  }
}

template<typename It, typename Compare =  std::less<>>
void MergeSortInPlace(It begin, It end, Compare compare = Compare{}) {
  if(begin < std::prev(end, 1)) {
    auto mid = std::next(begin, std::distance(begin, end)/2);
    MergeSortInPlace(begin, mid, compare);
    MergeSortInPlace(mid, end, compare);
    MergeSortInPlaceUtil(begin, mid, end, compare);
  }
}

template<typename It, typename Compare = std::less<>>
It QuickSortPartitionUtil(It begin, It last, Compare compare = Compare{}) {
  auto pivot = last;
  auto current = begin;
  for(auto it = begin; it < last; ++it) {
    if(compare(*it, *pivot)) {
      std::iter_swap(current, it);
      ++current;
    }
  }
  std::iter_swap(current, pivot);
  return current;
}

template<typename It, typename Compare = std::less<>>
void QuickSortUtil(It begin, It last, Compare compare = Compare{}) {
  if(begin < last) {
    auto partition = QuickSortPartitionUtil(begin, last, compare);
    QuickSortUtil(begin, std::prev(partition, 1), compare);
    QuickSortUtil(std::next(partition, 1), last, compare);
  }
}

template<typename It, typename Compare = std::less<>>
void QuickSort(It begin, It end, Compare compare = Compare{}) {
    QuickSortUtil(begin, std::prev(end), compare);
}

template<typename It, typename Compare = std::less<>>
void HeapSortHeapifyUtil(It begin, It current, It end, Compare compare) {
  auto it = current;
  auto dist = std::distance(begin, current)*2;
  auto l = std::next(begin, dist + 1);
  auto r = std::next(l, 1);
  if(l < end && *l > *it) it = l;
  if(r < end && *r > *it) it = r;
  if(it != current) {
    std::iter_swap(it, current);
    HeapSortHeapifyUtil(begin, it, end, compare);
  }
}

template<typename It, typename Compare = std::less<>>
void HeapSort(It begin, It end, Compare compare = Compare{}) {
  for(auto it = std::next(begin, std::distance(begin, end)/2 - 1); it != std::prev(begin); --it) {
    HeapSortHeapifyUtil(begin, it, end, compare);
  }
  for(auto it = std::prev(end); it != begin; --it) {
    std::iter_swap(begin, it);
    HeapSortHeapifyUtil(begin, begin, it, compare);
  }
}

}  // namespace sort
