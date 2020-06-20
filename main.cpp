#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "binary_tree.hpp"


std::shared_ptr<binary_tree::Node> func_bt_1() {
  std::uniform_int_distribution<> dist(5, 10);
  return binary_tree::generateRandomBinaryTree(5, dist);
}

std::shared_ptr<binary_tree::Node> func_bt_2() {
  std::uniform_int_distribution<> dist(0, 4);
  return binary_tree::generateRandomBinaryTree(5, dist);
}

int main() {

  // Test inOrder, preOrder(), postOrder()
  std::shared_ptr<binary_tree::Node> root = std::make_shared<binary_tree::Node>(0);
  root->left_ = std::make_shared<binary_tree::Node>(1);
  root->right_ = std::make_shared<binary_tree::Node>(2);
  root->left_->left_ = std::make_shared<binary_tree::Node>(3);
  root->left_->right_ = std::make_shared<binary_tree::Node>(4);

  std::cout << '\n' << "Pre order: " << '\n';
  binary_tree::preOrder(root);

  std::cout << '\n' << "In order: " << '\n';
  binary_tree::inOrder(root);

  std::cout << '\n' << "Post order: " << '\n';
  binary_tree::postOrder(root);
  std::cout << '\n';

  // Test recondtruction using constructBinaryTreeFromPreIn() and constructBinaryTreeFromInPost()
  std::vector<double> pre = {0, 1, 3, 4, 2};
  std::vector<double> in = {3, 1, 4, 0, 2};
  std::vector<double> post = {3, 4, 1, 2, 0};

  auto node_pre_in = binary_tree::constructBinaryTreeFromPreIn(pre, in, 0, in.size()-1);
  std::cout << '\n' << "Reconstructed binary tree pre order: " << '\n';
  binary_tree::preOrder(node_pre_in);
  std::cout << '\n';

  auto node_in_post = binary_tree::constructBinaryTreeFromInPost(in, post, 0, in.size() - 1);
  std::cout << '\n' << "Reconstructed binary tree pre order: " << '\n';
  binary_tree::preOrder(node_in_post);
  std::cout << '\n';

  // Test generation of binary tree using generateRandomBinaryTree(depth, min, max)
  auto node_random_bt = binary_tree::generateRandomBinaryTree(5, 0, 1);
  std::cout << '\n' << "[generateRandomBinaryTree(depth, min, max)] Generated random binary tree pre order: " << '\n';
  binary_tree::preOrder(node_random_bt);
  std::cout << '\n';

  // Test generation of full binary tree using generateRandomFullBinaryTree(depth, min, max)
  auto node_random_full_bt = binary_tree::generateRandomFullBinaryTree(5, 2, 3);
  std::cout << '\n' << "[generateRandomFullBinaryTree(depth, min, max)] Generated random full binary tree pre order: " << '\n';
  binary_tree::preOrder(node_random_full_bt);
  std::cout << '\n';

  // Test generation of full binary tree using generateRandomFullBinaryTree(depth, distribution)
  std::uniform_int_distribution<> d(5, 10);
  auto node_random_dist_bt = binary_tree::generateRandomBinaryTree(5, d);
  std::cout << '\n' << "[generateRandomFullBinaryTree(depth, distribution)] Generated random binary tree pre order: " << '\n';
  binary_tree::preOrder(node_random_dist_bt);
  std::cout << '\n';

  // Test generation of full binary tree using generateRandomFullBinaryTree(depth, distribution)
  std::uniform_real_distribution<> df(5, 10);
  auto node_random_dist_full_bt = binary_tree::generateRandomBinaryTree(5, df);
  std::cout << '\n' << "[generateRandomFullBinaryTree(depth, distribution)] Generated random full binary tree pre order: " << '\n';
  binary_tree::preOrder(node_random_dist_full_bt);
  std::cout << '\n';

  // Test generation of binary trees when being created simultaneously
  // There were sme static variables used previously that would obviously cause issues in multithreading
  // They have since been removed
  std::future<std::shared_ptr<binary_tree::Node>> bt_1 = std::async(std::launch::async, func_bt_1);
  std::future<std::shared_ptr<binary_tree::Node>> bt_2 = std::async(std::launch::async, func_bt_2);
  std::cout << '\n' << "[multithreading] Generated random full binary trees pre order: " << '\n';
  binary_tree::preOrder(bt_1.get());
  std::cout << '\n';
  binary_tree::preOrder(bt_2.get());
  std::cout << '\n';

  return 0;
}
