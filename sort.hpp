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



//int main() {
//	std::vector<int> to_sort{10,8,9,3,5,6,1,4,2,7};
//	BubbleSort(to_sort);
//	InsertionSort(to_sort);
//	InsertionSortWithBinarySearch(to_sort);
//	for(const auto& ele : to_sort) {
//		std::cout << ele << ' ';
//	}
//	std::cout << '\n';
//	return 0;
//}
