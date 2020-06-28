#include <iostream>
#include <stack>

#include "binary_search_tree.hpp"

namespace binary_search_tree {

std::shared_ptr<binary_tree::Node> constructTreeFromPreOrderUtil_N2(
  const std::vector<int>& pre_order, int start, int end, int& pre_index) {

  const int s = pre_order.size();
  if(start > end || pre_index >= s) {
    return nullptr;
  }

  std::shared_ptr<binary_tree::Node> node =
    std::make_shared<binary_tree::Node>(pre_order[pre_index]);
  pre_index++;

  if(start != end) {
    int break_index;
    for(int i = start; i <= end; i++) {
      if(pre_order[i] > node->value_) {
        break_index = i;
        break;
      }
    }
    node->left_ = constructTreeFromPreOrderUtil_N2(pre_order, pre_index, break_index - 1, pre_index);
    node->right_ = constructTreeFromPreOrderUtil_N2(pre_order, break_index, end, pre_index);
  }
  return node;
}

std::shared_ptr<binary_tree::Node> constructTreeFromPreOrderUtil_N(
  const std::vector<double>& pre_order, double min, double max, int& pre_index) {

  const int s = pre_order.size();
  if(pre_index > s) {
    return nullptr;
  }

  std::shared_ptr<binary_tree::Node> node =  nullptr;
  if(pre_order[pre_index] > min && pre_order[pre_index] < max) {
    node = std::make_shared<binary_tree::Node>(pre_order[pre_index]);
    pre_index++;

    if(pre_index < s) {
      node->left_ = constructTreeFromPreOrderUtil_N(pre_order, min, node->value_, pre_index);
      node->right_ = constructTreeFromPreOrderUtil_N(pre_order, node->value_, max, pre_index);
    }
  }

  return node;
}

// max is the value of the node above
std::shared_ptr<binary_tree::Node> constructTreeFromPreOrderUtilStack_N(const std::vector<double>& pre_order) {
  int s = pre_order.size();
  if(s==0) return nullptr;

  int pre_index = 0;
  std::stack<std::shared_ptr<binary_tree::Node>> node_stack;
  std::cout << "main creating tree node with value " << pre_order[pre_index] << '\n';
  std::shared_ptr<binary_tree::Node> root = std::make_shared<binary_tree::Node>(pre_order[pre_index]);
  node_stack.push(root);
  pre_index++;
  while(pre_index < s) {
    std::shared_ptr<binary_tree::Node> node = nullptr;
    while(!node_stack.empty() && node_stack.top()->value_ <= pre_order[pre_index]) {
      node = node_stack.top();
      node_stack.pop();
    }
    if(node!=nullptr) {
      std::cout << "right creating tree node with value " << pre_order[pre_index] << '\n';
      node->right_ = std::make_shared<binary_tree::Node>(pre_order[pre_index]);
      node_stack.push(node->right_);
    } else {
      node = node_stack.top();
      std::cout << "left creating tree node with value " << pre_order[pre_index] << '\n';
      node->left_ = std::make_shared<binary_tree::Node>(pre_order[pre_index]);
      node_stack.push(node->left_);
    }
    pre_index++;
  }
  return root;
}

std::shared_ptr<binary_tree::Node> constructTreeFromPreOrder(const std::vector<double>& pre_order) {
  const int s = pre_order.size();
  int pre_index = 0;
  return constructTreeFromPreOrderUtilStack_N(pre_order);
}

std::shared_ptr<binary_tree::Node> findLowestCommonAncestorUtil(const std::shared_ptr<binary_tree::Node>& node, double v1, double v2) {
  if (node == nullptr) return nullptr;
  else if (node->value_ > v1 && node->value_ < v2 ) return node;
  else if (node->value_ > v2) return findLowestCommonAncestorUtil(node->left_, v1, v2);
  else if (node->value_ < v1) return findLowestCommonAncestorUtil(node->right_, v1, v2);
  else return nullptr;
}

std::shared_ptr<binary_tree::Node> findLowestCommonAncestor(const std::shared_ptr<binary_tree::Node>& root, double v1, double v2) {
  if(v2 < v1) std::swap(v1, v2);
  std::shared_ptr<binary_tree::Node> node = findLowestCommonAncestorUtil(root, v1, v2);
  if(node == nullptr) {
    std::cout << __FUNCTION__ << " | " << "No ancestor found" << '\n';
  } else {
    std::cout << __FUNCTION__ << " | " << "Lowest common ancestor has value: " << node->value_ << '\n';
  }
  return node;
}

}; // namespace binary_search_tree

//
//
// int main() {
//   std::vector<double> pre = {10, 5, 1, 7, 40, 50};
//   std::shared_ptr<binary_tree::Node> root = binary_search_tree::constructTreeFromPreOrder(pre);
//   std::cout << "The reconstructed tree is: " << '\n';
//   inOrder(root);
//   std::cout << '\n';
//   return 0;
// }
