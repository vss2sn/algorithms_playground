#ifndef BINARY_TREE_HPP
#define BINARY_TREE_HPP

#include <memory>
#include <vector>
#include <random>

namespace binary_tree {

class Node {
public:
  Node(double value = 0);
  std::shared_ptr<Node> left_ = nullptr, right_ = nullptr;
  double value_ = 0;
};

void inOrder(const std::shared_ptr<Node>& node);
void postOrder(const std::shared_ptr<Node>& node);
void preOrder(const std::shared_ptr<Node>& node);

std::shared_ptr<Node> constructBinaryTreeFromPreIn(const std::vector<double>& pre_trav, const std::vector<double>& in_trav, const int start, const int end);
std::shared_ptr<Node> constructBinaryTreeFromInPost(const std::vector<double>& in_trav, const std::vector<double>& post_trav, const int start, const int end);

std::shared_ptr<Node> constructBinaryTreeFromInPostUtil(const std::vector<double>& in_trav, const std::vector<double>& post_trav, const int start, const int end, int& post_index);
std::shared_ptr<Node> constructBinaryTreeFromPreInUtil(const std::vector<double>& pre_trav, const std::vector<double>& in_trav, const int start, const int end, int& pre_index);

std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, const double& min_val, const double& max_val);
std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, const double& min_val, const double& max_val);

template <typename T>
std::shared_ptr<Node> generateRandomBinaryTreeUtil(const int max_depth, T& dist, std::mt19937& gen, double& random_null_thresh);

template<typename T>
std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, T& dist);

template <typename T>
std::shared_ptr<Node> generateRandomFullBinaryTreeUtil(const int max_depth, T& dist, std::mt19937& gen, double& random_null_thresh);

template<typename T>
std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, T& dist);

// Templated function definitions:

template <typename T>
std::shared_ptr<Node> generateRandomFullBinaryTreeUtil(const int max_depth, T& dist, std::mt19937& gen, double& random_null_thresh) {

  // Add check to make sure that the distribution passed in makes sense? or let SFINAE deal with that?
  if (max_depth == 0) return nullptr;
  double value = dist(gen);
  std::shared_ptr<Node> node = std::make_shared<Node>(value);

  if (random_null_thresh < static_cast<double>(std::rand())/RAND_MAX) {
    node->left_ = generateRandomFullBinaryTreeUtil(max_depth - 1, dist, gen, random_null_thresh);
    node->right_ = generateRandomFullBinaryTreeUtil(max_depth - 1, dist, gen, random_null_thresh);
  } else {  // This is not really required but attempting to future proof
    node->left_ = nullptr;
    node->right_ = nullptr;
  }
  return node;
}

template<typename T>
std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, T& dist) {
  if (max_depth == 0) return nullptr;
  std::random_device rd;
  std::mt19937 gen(rd());

  std::shared_ptr<Node> node;
  double random_null_thresh = 1.0/max_depth;
  node = generateRandomFullBinaryTreeUtil(max_depth, dist, gen, random_null_thresh);

  return node;
}

template <typename T>
std::shared_ptr<Node> generateRandomBinaryTreeUtil(const int max_depth, T& dist, std::mt19937& gen, double& random_null_thresh) {

  // Add check to make sure that the distribution passed in makes sense? or let SFINAE deal with that?
  if (max_depth == 0) return nullptr;
  double value = dist(gen);
  std::shared_ptr<Node> node = std::make_shared<Node>(value);

  if (random_null_thresh < static_cast<double>(std::rand())/RAND_MAX) {
    node->left_ = generateRandomBinaryTreeUtil(max_depth - 1, dist, gen, random_null_thresh);
  } else {  // This is not really required but attempting to future proof
    node->left_ = nullptr;
  }

  if (random_null_thresh < static_cast<double>(std::rand())/RAND_MAX) {
    node->right_ = generateRandomBinaryTreeUtil(max_depth - 1, dist, gen, random_null_thresh);
  } else {  // This is not really required but attempting to future proof
    node->right_ = nullptr;
  }

  return node;
}

template<typename T>
std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, T& dist) {
  if (max_depth == 0) return nullptr;

  std::random_device rd;
  std::mt19937 gen(rd());

  double random_null_thresh = 1.0/max_depth;

  return generateRandomBinaryTreeUtil(max_depth, dist, gen, random_null_thresh);
}

}; // namespace binary_tree

#endif  // BINARY_TREE_HPP
