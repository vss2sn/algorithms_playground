#include <algorithm>
#include <iostream>
#include <vector>

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
		// use std copy backward?
		// while(j >=0 && to_sort[j] > key) {
		// 	j--;
		// }
		// std::copy_backward(&to_sort[j+1], &to_sort[i], &to_sort[i+1]);
		while (j >= 0 && to_sort[j] > key) {
			to_sort[j+1] = to_sort[j];
			j--;
		}
		to_sort[j+1] = key;
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

		// use std copy backaward
		//std::copy_backward(&to_sort[loc], &to_sort[i], &to_sort[i+1]);
		while(j >= loc) {
			to_sort[j+1] = to_sort[j];
			j--;
		}
		to_sort[loc] = key;
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
