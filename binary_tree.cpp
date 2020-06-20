#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "binary_tree.hpp"

namespace binary_tree {

Node::Node(double value) {
  value_ = value;
}

void inOrder(const std::shared_ptr<Node>& node) {
  if (node == nullptr) return;
  inOrder(node->left_);
  std::cout << node->value_ << ' ';
  inOrder(node->right_);
}

void postOrder(const std::shared_ptr<Node>& node) {
  if (node == nullptr) return;
  postOrder(node->left_);
  postOrder(node->right_);
  std::cout << node->value_ << ' ';
}

void preOrder(const std::shared_ptr<Node>& node) {
  if (node == nullptr) return;
  std::cout << node->value_ << ' ';
  preOrder(node->left_);
  preOrder(node->right_);
}


// Not using std::set for the inorder traversal to allow non-unique elements
std::shared_ptr<Node> constructBinaryTreeFromPreIn(
  const std::vector<double>& pre_trav, const std::vector<double>& in_trav,
  const int start, const int end) {

  static int pre_index = 0;

  if (start > end) return nullptr;

  std::shared_ptr<Node> node = std::make_shared<Node>(pre_trav[pre_index]);
  int index = std::distance(in_trav.begin(), std::find(in_trav.begin() + start, in_trav.end() + end, node->value_));

  pre_index++;
  if (start == end) return node;

  node->left_ = constructBinaryTreeFromPreIn(pre_trav, in_trav, start, index - 1);
  node->right_ = constructBinaryTreeFromPreIn(pre_trav, in_trav, index + 1, end);

  return node;
}

// Not using std::set for the inorder traversal to allow non-unique elements
std::shared_ptr<Node> constructBinaryTreeFromInPost(
  const std::vector<double>& in_trav, const std::vector<double>& post_trav,
  const int start, const int end) {

  static int post_index = in_trav.size()-1;

  if (start > end) return nullptr;

  std::shared_ptr<Node> node = std::make_shared<Node>(post_trav[post_index]);
  int index = std::distance(in_trav.begin(), std::find(in_trav.begin() + start, in_trav.end() + end + 1, node->value_));

  post_index--;
  if (start == end) return node;

  node->right_ = constructBinaryTreeFromInPost(in_trav, post_trav, index + 1, end);
  node->left_ = constructBinaryTreeFromInPost(in_trav, post_trav, start, index - 1);

  return node;
}

std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, const double min_val, const double max_val) {
  if(max_depth == 0) return nullptr;
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<> dist(min_val, max_val);
  static double random_null = dist(gen);
  double value = dist(gen);
  std::shared_ptr<Node> node = std::make_shared<Node>(value);
  if(random_null < dist(gen)) {
    node->left_ = generateRandomBinaryTree(max_depth - 1, min_val, max_val);
  } else {
    node->left_ = nullptr;
  }
  if(random_null < dist(gen)) {
    node->right_ = generateRandomBinaryTree(max_depth - 1, min_val, max_val);
  } else {
    node->right_ = nullptr;
  }
  return node;
}

std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, const double min_val, const double max_val) {
  if(max_depth == 0) return nullptr;
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<> dist(min_val, max_val);
  static double random_null = dist(gen);
  double value = dist(gen);
  std::shared_ptr<Node> node = std::make_shared<Node>(value);

  if(random_null < dist(gen)) {
    node->left_ = generateRandomFullBinaryTree(max_depth - 1, min_val, max_val);
    node->right_ = generateRandomFullBinaryTree(max_depth - 1, min_val, max_val);
  } else {
    node->left_ = nullptr;
    node->right_ = nullptr;
  }
  return node;
}
//
// std::shared_ptr<Node> generateRandomCompleteBinaryTree(const int max_depth, const double min_val, const double max_val) {
//   if(max_depth == 0) return nullptr;
//   static std::random_device rd;
//   static std::mt19937 gen(rd());
//   static std::uniform_real_distribution<> dist(min_val, max_val);
//   static double random_null = dist(gen);
//   double value = dist(gen);
//   std::shared_ptr<Node> node = std::make_shared<Node>(value);
//
//   if(random_null < dist(gen)) {
//     node->left_ = generateRandomBinaryTree(max_depth - 1, min_val, max_val);
//     node->right_ = generateRandomBinaryTree(max_depth - 1, min_val, max_val);
//   } else {
//     node->left_ = nullptr;
//     node->right_ = nullptr;
//   }
//   return node;
// }
}; // namespace binary_tree
