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
  MergeSort(begin, mid);
  MergeSort(mid, end);
  MergeSortMergeUtil(begin, mid, end, compare);
}

}  // namespace sort
