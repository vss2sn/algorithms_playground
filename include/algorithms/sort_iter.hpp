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
    std::iter_swap(it, std::min_element(it, end));
  }
}

}  // namespace sort
