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

// std::shared_ptr<binary_tree::Node> func_bt_1() {
//   std::uniform_int_distribution<> dist(5, 10);
//   return binary_tree::generateRandomBinaryTree(5, dist);
// }
//
// std::shared_ptr<binary_tree::Node> func_bt_2() {
//   std::uniform_int_distribution<> dist(0, 4);
//   return binary_tree::generateRandomBinaryTree(5, dist);
// }
// 
// int main() {
//
//   // Test inOrder, preOrder(), postOrder()
//   std::shared_ptr<binary_tree::Node> root = std::make_shared<binary_tree::Node>(0);
//   root->left_ = std::make_shared<binary_tree::Node>(1);
//   root->right_ = std::make_shared<binary_tree::Node>(2);
//   root->left_->left_ = std::make_shared<binary_tree::Node>(3);
//   root->left_->right_ = std::make_shared<binary_tree::Node>(4);
//
//   std::cout << '\n' << "Pre order: " << '\n';
//   binary_tree::preOrder(root);
//
//   std::cout << '\n' << "In order: " << '\n';
//   binary_tree::inOrder(root);
//
//   std::cout << '\n' << "Post order: " << '\n';
//   binary_tree::postOrder(root);
//   std::cout << '\n';
//
//   // Test recondtruction using constructBinaryTreeFromPreIn() and constructBinaryTreeFromInPost()
//   std::vector<double> pre = {0, 1, 3, 4, 2};
//   std::vector<double> in = {3, 1, 4, 0, 2};
//   std::vector<double> post = {3, 4, 1, 2, 0};
//
//   auto node_pre_in = binary_tree::constructBinaryTreeFromPreIn(pre, in, 0, in.size()-1);
//   std::cout << '\n' << "Reconstructed binary tree pre order: " << '\n';
//   binary_tree::preOrder(node_pre_in);
//   std::cout << '\n';
//
//   auto node_in_post = binary_tree::constructBinaryTreeFromInPost(in, post, 0, in.size() - 1);
//   std::cout << '\n' << "Reconstructed binary tree pre order: " << '\n';
//   binary_tree::preOrder(node_in_post);
//   std::cout << '\n';
//
//   // Test generation of binary tree using generateRandomBinaryTree(depth, min, max)
//   auto node_random_bt = binary_tree::generateRandomBinaryTree(5, 0, 1);
//   std::cout << '\n' << "[generateRandomBinaryTree(depth, min, max)] Generated random binary tree pre order: " << '\n';
//   binary_tree::preOrder(node_random_bt);
//   std::cout << '\n';
//
//   // Test generation of full binary tree using generateRandomFullBinaryTree(depth, min, max)
//   auto node_random_full_bt = binary_tree::generateRandomFullBinaryTree(5, 2, 3);
//   std::cout << '\n' << "[generateRandomFullBinaryTree(depth, min, max)] Generated random full binary tree pre order: " << '\n';
//   binary_tree::preOrder(node_random_full_bt);
//   std::cout << '\n';
//
//   // Test generation of full binary tree using generateRandomFullBinaryTree(depth, distribution)
//   std::uniform_int_distribution<> d(5, 10);
//   auto node_random_dist_bt = binary_tree::generateRandomBinaryTree(5, d);
//   std::cout << '\n' << "[generateRandomFullBinaryTree(depth, distribution)] Generated random binary tree pre order: " << '\n';
//   binary_tree::preOrder(node_random_dist_bt);
//   std::cout << '\n';
//
//   // Test generation of full binary tree using generateRandomFullBinaryTree(depth, distribution)
//   std::uniform_real_distribution<> df(5, 10);
//   auto node_random_dist_full_bt = binary_tree::generateRandomBinaryTree(5, df);
//   std::cout << '\n' << "[generateRandomFullBinaryTree(depth, distribution)] Generated random full binary tree pre order: " << '\n';
//   binary_tree::preOrder(node_random_dist_full_bt);
//   std::cout << '\n';
//
//   // Test generation of binary trees when being created simultaneously
//   // There were sme static variables used previously that would obviously cause issues in multithreading
//   // They have since been removed
//   std::future<std::shared_ptr<binary_tree::Node>> bt_1 = std::async(std::launch::async, func_bt_1);
//   std::future<std::shared_ptr<binary_tree::Node>> bt_2 = std::async(std::launch::async, func_bt_2);
//   std::cout << '\n' << "[multithreading] Generated random full binary trees pre order: " << '\n';
//   binary_tree::preOrder(bt_1.get());
//   std::cout << '\n';
//   binary_tree::preOrder(bt_2.get());
//   std::cout << '\n';
//
//   return 0;
// }
