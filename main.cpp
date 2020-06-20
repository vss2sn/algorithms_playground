#include <future>
#include <iostream>
#include <memory>
#include <thread>

#include "binary_tree.hpp"


std::shared_ptr<binary_tree::Node> func_bt_1() {
  return binary_tree::generateRandomBinaryTree(30,0,1);
}

std::shared_ptr<binary_tree::Node> func_bt_2() {
  return binary_tree::generateRandomBinaryTree(30,2,3);
}

int main() {

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

  auto node_random_bt = binary_tree::generateRandomBinaryTree(5, 0, 1);
  std::cout << '\n' << "Generated random binary tree pre order: " << '\n';
  binary_tree::preOrder(node_random_bt);
  std::cout << '\n';

  auto node_random_full_bt = binary_tree::generateRandomFullBinaryTree(5, 2, 3);
  std::cout << '\n' << "Generated random full binary tree pre order: " << '\n';
  binary_tree::preOrder(node_random_full_bt);
  std::cout << '\n';

  std::future<std::shared_ptr<binary_tree::Node>> bt_1 = std::async(std::launch::async, func_bt_1);
  std::future<std::shared_ptr<binary_tree::Node>> bt_2 = std::async(std::launch::async, func_bt_2);
  std::cout << '\n';
  binary_tree::preOrder(bt_1.get());
  std::cout << '\n';
  std::cout << '\n';
  binary_tree::preOrder(bt_2.get());
  std::cout << '\n';
  return 0;
}
