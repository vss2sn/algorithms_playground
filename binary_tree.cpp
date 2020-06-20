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
std::shared_ptr<Node> constructBinaryTreeFromPreIn(const std::vector<double>& pre_trav, const std::vector<double>& in_trav,
const int start, const int end) {
  int pre_index = 0;
  return constructBinaryTreeFromPreInUtil(pre_trav, in_trav, start, end, pre_index);
}

std::shared_ptr<Node> constructBinaryTreeFromPreInUtil(
  const std::vector<double>& pre_trav, const std::vector<double>& in_trav,
  const int start, const int end, int& pre_index) {

  if (start > end) return nullptr;

  std::shared_ptr<Node> node = std::make_shared<Node>(pre_trav[pre_index]);
  int index = std::distance(in_trav.begin(), std::find(in_trav.begin() + start, in_trav.end() + end, node->value_));

  pre_index++;
  if (start == end) return node;

  node->left_ = constructBinaryTreeFromPreInUtil(pre_trav, in_trav, start, index - 1, pre_index);
  node->right_ = constructBinaryTreeFromPreInUtil(pre_trav, in_trav, index + 1, end, pre_index);

  return node;
}

// Not using std::set for the inorder traversal to allow non-unique elements
std::shared_ptr<Node> constructBinaryTreeFromInPost(
  const std::vector<double>& in_trav, const std::vector<double>& post_trav,
  const int start, const int end) {
    int post_index = in_trav.size()-1;
    return constructBinaryTreeFromInPostUtil(in_trav, post_trav, start, end, post_index);
}

std::shared_ptr<Node> constructBinaryTreeFromInPostUtil(
  const std::vector<double>& in_trav, const std::vector<double>& post_trav,
  const int start, const int end, int& post_index) {

  if (start > end) return nullptr;

  std::shared_ptr<Node> node = std::make_shared<Node>(post_trav[post_index]);
  int index = std::distance(in_trav.begin(), std::find(in_trav.begin() + start, in_trav.end() + end + 1, node->value_));

  post_index--;
  if (start == end) return node;

  node->right_ = constructBinaryTreeFromInPostUtil(in_trav, post_trav, index + 1, end, post_index);
  node->left_ = constructBinaryTreeFromInPostUtil(in_trav, post_trav, start, index - 1, post_index);

  return node;
}

std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, const double& min_val, const double& max_val) {
  if (max_depth == 0) return nullptr;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(min_val, max_val);
  double random_null_thresh = 1.0/max_depth;
  return generateRandomBinaryTreeUtil(max_depth, dist, gen, random_null_thresh);
}

std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, const double& min_val, const double& max_val) {
  if (max_depth == 0) return nullptr;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(min_val, max_val);
  double random_null_thresh = 1.0/max_depth;
  return generateRandomFullBinaryTreeUtil(max_depth, dist, gen, random_null_thresh);
}

}; // namespace binary_tree
