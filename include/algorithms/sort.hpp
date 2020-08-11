#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace sort {

template<typename T>
void BubbleSort(std::vector<T>& to_sort) {
	int v = to_sort.size();

	for(int i = 0; i < v; i++) {
		bool swapped = false;
		for(int j=0; j < v - i - 1; j++) {
			if(to_sort[j] > to_sort[j+1]) {
				std::swap(to_sort[j], to_sort[j+1]);
				swapped = true;
			}
		}
		if(!swapped) break;
	}
}

template<typename T>
void InsertionSort(std::vector<T>& to_sort) {
	int v = to_sort.size();

	for (int i=1; i< v; i++) {
		T key = to_sort[i];
		int j = i-1;
		while (j >= 0 && to_sort[j] > key) {
			to_sort[j+1] = to_sort[j];
			j--;
		}
		to_sort[j+1] = key;
	}
}

// Playing with STL
template<typename T>
void InsertionSortWithCb(std::vector<T>& to_sort) {
	int v = to_sort.size();

	for (int i=1; i< v; i++) {
		T key = to_sort[i];
		int j = i-1;
		while(j >=0 && to_sort[j] > key) {
			j--;
		}
		std::copy_backward(&to_sort[j+1], &to_sort[i], &to_sort[i+1]);
		to_sort[j+1] = key;
	}
}

// Playing with STL
template<typename T>
void InsertionSortWithRiCb(std::vector<T>& to_sort) {
	int v = to_sort.size();

	for (int i=1; i< v; i++) {
		T key = to_sort[i];
		int j = std::distance(
			std::find_if(to_sort.rbegin() + v-i, to_sort.rend(), [&](T temp){return temp <= key;}),
			to_sort.rend()
		);
		if(j+1 <= i) {
			std::copy_backward(&to_sort[j], &to_sort[i], &to_sort[i+1]);
			to_sort[j] = key;
		}
	}
}


template<typename T>
int BinarySearchUtil(const std::vector<T>& to_search, T key, int begin, int end) {
	if(begin >= end) {
		if(key > to_search[begin]) {
			return begin + 1;
		} else {
			return begin;
		}
	}

	int mid = (begin + end)/2;
	if(key == to_search[mid]) {
		return	mid + 1;
	} else if(key > to_search[mid]) {
		return BinarySearchUtil(to_search, key, mid + 1, end);
	} else {
		return BinarySearchUtil(to_search, key, begin, mid - 1);
	}
}

template<typename T>
void InsertionSortWithBinarySearch(std::vector<T>& to_sort) {
	int v = to_sort.size();

	for(int i = 1; i < v; i++) {
		int j = i-1;
		T key = to_sort[i];
		int loc = BinarySearchUtil(to_sort, key, 0, j);

		while(j >= loc) {
			to_sort[j+1] = to_sort[j];
			j--;
		}
		to_sort[loc] = key;
	}
}

// Playing with STL
// NOTE: This is faster than the `STL optimized version` below
template<typename T>
void InsertionSortWithBinarySearchWithCb(std::vector<T>& to_sort) {
	int v = to_sort.size();

	for(int i = 1; i < v; i++) {
		int j = i-1;
		T key = to_sort[i];
		int loc = BinarySearchUtil(to_sort, key, 0, j);
		std::copy_backward(&to_sort[loc], &to_sort[i], &to_sort[i+1]);
		to_sort[loc] = key;
	}
}

// Playing with STL
template<typename T>
void InsertionSortOptimized(std::vector<T>& to_sort) {
	for(auto it = to_sort.begin(); it != to_sort.end(); ++it) {
		std::rotate(std::upper_bound(to_sort.begin(), it, *it), it, it + 1);
	}
}

// NOTE: This is faster than the `STL optimized version` above
template<typename T>
void InsertionSortOptimizedWithCb(std::vector<T>& to_sort) {
	for(auto it = to_sort.begin() + 1; it != to_sort.end(); ++it) {
		const T key = *it;
		const auto it_next_larger = std::upper_bound(to_sort.begin(), it, *it);
		std::copy_backward(it_next_larger, it, it+1);
		*it_next_larger = key;
	}
}

template<typename T>
void MergeSortMergeUtil(std::vector<T>& to_sort, int begin, int mid, int end) {
	int i = 0, j = 0, k = begin;

	// NOTE: the constructor is vector v[a, b); b not inclusive
	std::vector<T> left(to_sort.begin() + begin, to_sort.begin() + mid + 1);
	std::vector<T> right(to_sort.begin() + mid + 1, to_sort.begin() + end + 1);

	while(i < left.size() && j < right.size()) {
		if(left[i] < right[j]) {
			to_sort[k] = left[i];
			i++;
		} else {
			to_sort[k] = right[j];
			j++;
		}
		k++;
	}
	while(j < right.size()) {
		to_sort[k] = right[j];
		j++;
		k++;
	}
	while(i < left.size()) {
		to_sort[k] = left[i];
		i++;
		k++;
	}
}

template<typename T>
void MergeSortUtil(std::vector<T>& to_sort, int begin, int end) {
	if(begin < end) {
		int mid = begin + (end - begin)/2;
		MergeSortUtil(to_sort, begin, mid);
		MergeSortUtil(to_sort, mid+1, end);
		MergeSortMergeUtil(to_sort, begin, mid, end);
	}
}

template<typename T>
void MergeSort(std::vector<T>& to_sort) {
	int v = to_sort.size();
	MergeSortUtil(to_sort, 0, v-1);
}

template<typename T>
void MergeSortMergeInPlaceUtil(std::vector<T>& to_sort, int begin, int mid, int end) {
	int begin2 = mid+1;

	while(begin <= mid && begin2 <= end) {
		if(to_sort[begin] <= to_sort[begin2]) {
			++begin;
		} else {
			// std::rotate(to_sort.begin() + begin,
			// 						to_sort.begin() + begin2,
			// 						to_sort.begin() + begin2 + 1);
			// Alternate syntax as this file rarely uses iterators
			std::rotate(&to_sort[begin], &to_sort[begin2], &to_sort[begin2+1]);
			// NOTE: The commented out section below is also correct
			// and is the cb version which appears to be faster
			// const auto key = *(to_sort.begin() + begin2);
			// std::copy_backward(to_sort.begin() + begin,
			// 									 to_sort.begin() + begin2,
			// 									 to_sort.begin() + begin2+1);
		  // *(to_sort.begin() + begin) = key;
			++begin;
			++mid;
			++begin2;
		}
	}
}

template<typename T>
void MergeSortInPlaceUtil(std::vector<T>& to_sort, int begin, int end) {
	if(begin < end) {
		int mid = begin + (end - begin)/2;
		MergeSortUtil(to_sort, begin, mid);
		MergeSortUtil(to_sort, mid+1, end);
		MergeSortMergeInPlaceUtil(to_sort, begin, mid, end);
	}
}

template<typename T>
void MergeSortInPlace(std::vector<T>& to_sort) {
	int v = to_sort.size();
	MergeSortInPlaceUtil(to_sort, 0, v-1);
}

template<typename T>
int QuickSortPartition(std::vector<T>& to_sort, int start, int end) {
	T pivot = to_sort[end];
	int i = start;
	for(int j = start; j < end; j++) {
		if(to_sort[j] < pivot) {
			std::swap(to_sort[i], to_sort[j]);
			i++;
		}
	}
	std::swap(to_sort[end], to_sort[i]);
	return i;
}

template<typename T>
void QuickSortUtil(std::vector<T>& to_sort, int start, int end) {
	if(start < end) {
		int p = QuickSortPartition(to_sort, start, end);
		QuickSortUtil(to_sort, start, p-1);
		QuickSortUtil(to_sort, p+1, end);
	}
}

template<typename T>
void QuickSort(std::vector<T>& to_sort) {
	int v = to_sort.size();
	QuickSortUtil(to_sort, 0, v-1);
}

template<typename T>
void SelectionSort(std::vector<T>& to_sort) {
  int v = to_sort.size();
  for(int i = 0; i < v - 1; i++) {
    int min_index = i;
    for(int j = i+1; j < v; j++) {
      if(to_sort[min_index] > to_sort[j]) {
        min_index = j;
      }
    }
		std::swap(to_sort[min_index], to_sort[i]);
  }
}

// Only valid for positive values
template <typename T>
void BucketSort(std::vector<T>& to_sort, const int n = 10) {
	if(n <= 0) {
		std::cout << "Please set n >= 1" << '\n';
		return;
	}
	int v = to_sort.size();
	std::vector<std::shared_ptr<std::vector<T>>> buckets(n);
	T key = *std::max_element(to_sort.begin(), to_sort.end());
	const double mp = (static_cast<double>(n)) / (key+1); // if key ==0 and key == n
	for(int i = 0; i<n; i++) {
		buckets[i] = std::make_shared<std::vector<T>>();
	}
	for(const auto& ele : to_sort) {
		buckets[std::floor(ele * mp)]->push_back(ele);
	}
	auto it = to_sort.begin();
	for(const auto& bucket : buckets) {
		QuickSort(*bucket);
		std::copy(bucket->begin(), bucket->end(), it);
		std::advance(it, bucket->size());
	}
}

template<typename T>
void CountingSort(std::vector<T>& to_sort) {
	const int v = to_sort.size();

	const T max = *std::max_element(to_sort.begin(), to_sort.end());
	const T min = *std::min_element(to_sort.begin(), to_sort.end());

	const int range = max - min + 1;

	std::vector<int> count(range, 0);
	for(const auto& ele:to_sort) {
		count[ele - min]++;
	}

	for(int i = 1; i < range; i++) {
		count[i]+=count[i-1];
	}

	std::vector<T> output(v);
	for(int i = v - 1; i >= 0; i--) {
		output[count[to_sort[i]-min]-1] = to_sort[i];
		count[to_sort[i]-min]--;
	}

	std::copy(output.begin(), output.end(), to_sort.begin());
}

template<typename T>
void HeapSortUtilHeapify(std::vector<T>& to_sort, const int n, const int i) {
	int largest = i;
	const int l = 2*i+1;
	const int r = 2*i+2;

	if(l < n && to_sort[l] > to_sort[largest]) {
		largest = l;
	}

	if(r < n && to_sort[r] > to_sort[largest]) {
		largest = r;
	}

	if(largest != i) {
		std::swap(to_sort[largest], to_sort[i]);
		HeapSortUtilHeapify(to_sort, n, largest);
	}
}

template<typename T>
void HeapSort(std::vector<T>& to_sort) {
	int v = to_sort.size();

	for(int i = v/2-1; i >= 0; i--){
		HeapSortUtilHeapify(to_sort, v, i);
	}

	for(int i = v-1; i >0; i--) {
		std::swap(to_sort[0], to_sort[i]);
		HeapSortUtilHeapify(to_sort, i, 0);
	}
}

}  // namespace sort
