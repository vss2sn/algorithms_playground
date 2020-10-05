# Algorithms playground #

#### This repository contains implementations of algorithms on graphs and binary trees, and sorting algorithms  ####

[![Build Status](https://travis-ci.com/vss2sn/algorithms_playground.svg?branch=master)](https://travis-ci.com/vss2sn/algorithms_playground)

#### To build and run: ####
     git clone https://github.com/vss2sn/algorithms_playground.git  
     cd algorithms_playground  
     mkdir build  
     cd build  
     cmake .. && make -j
     ./play

#### Algorithms: ####
1. Algorithms on graphs (using both Aadjacency lists and adjacency matrices)
   * Breadth first search
   * Depth first search
   * Prim's algorithms
   * Kruskal's algorithm
   * Ford Fulkerson algorithm
   * Find mother vertex
   * All paths between
   * Colour Graph (using the minimum number of colours)
   * Dijkstra's algorithm
   * Kosaraju's algorithm
   * Floyd Warshall algorithm
   * Bellman Ford algorithm
   * Find strongly connected components (using Kosaraju's algorithm)
2. Algorithms on binary trees
   * In order traversal
   * Pre order traversal
   * Post order traversal
   * construct binary tree from preorder and inorder traveral
   * construct binary tree from inorder and postorder traveral
3. Algorithms on binary search trees
   * construct BST from preorder traveral
   * insert into BST
   * find in BST
   * delete from BST
   * find lowest common ancestor
   * check pair sum BST
4. Sorting algorithms (using both indices and iterators)
   * Bubble sort
   * Insertion sort (5 variants including with binary search)
   * Merge sort
   * Merge sort in place
   * Quick sort
   * Quick sort in place
   * Selection sort
   * Bucket sort
   * Counting sort
   * Heap sort

#### TODO: ####
1. Add documentation for each algorithm
2. Add documentation explaining each algorithm
3. Add tests for each function
