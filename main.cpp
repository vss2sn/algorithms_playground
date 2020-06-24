#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "binary_tree.hpp"
#include "binary_search_tree.hpp"
#include "graphAM.hpp"

int main () {

  std::vector<double> pre = {10, 5, 1, 7, 40, 50};
  std::shared_ptr<binary_tree::Node> root = binary_search_tree::constructTreeFromPreOrder(pre);
  std::cout << "The reconstructed tree is: " << '\n';
  inOrder(root);
  std::cout << '\n';
  return 0;

}
