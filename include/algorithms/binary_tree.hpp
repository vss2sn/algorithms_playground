#ifndef BINARY_TREE_HPP
#define BINARY_TREE_HPP

#include <memory>
#include <vector>
#include <random>

namespace binary_tree {

/**
* @brief Class Node
* @details Represents a single tree node containing pointers to 2 leaf nodes and a value
*/
class Node {
public:
  Node(double value = 0);
  std::shared_ptr<Node> left_ = nullptr, right_ = nullptr;
  double value_ = 0;
};

/**
 * @brief Prints the inorder traversal of the tree
 * @param [in] node root node of tree/subtree
 */
void inOrder(const std::shared_ptr<Node>& node);

/**
 * @brief Prints the postorder traversal of the tree
 * @param [in] node root node of tree/subtree
 */
void postOrder(const std::shared_ptr<Node>& node);

/**
 * @brief Prints the preorder traversal of the tree
 * @param [in] node root node of tree/subtree
 */
void preOrder(const std::shared_ptr<Node>& node);

template<typename T>
void getInOrderUtil(const std::shared_ptr<Node>& node, T& in_order);

/**
 * @brief Prints the inorder traversal of the tree
 * @param [in] node root node of tree/subtree
 * @return templated data structure containing the iorder traversal of the tree
 */
template<typename T>
T getInOrder(const std::shared_ptr<Node>& node);

/**
 * @brief Construct a binary tree from preorder and inorder traversals
 * @param [in] pre_trav preorder traversal
 * @param [in] in_trav inorder traversal
 * @param [in] start start index of from which the tree is to be constructed
 * @param [in] end end index of from which the tree is to be constructed
 * @return pointer to root of newly constructed tree
 */
std::shared_ptr<Node> constructBinaryTreeFromPreIn(const std::vector<double>& pre_trav, const std::vector<double>& in_trav, const int start, const int end);

/**
 * @brief Construct a binary tree from inorder and postorder traversals
 * @param [in] in_trav inorder traversal
 * @param [in] post_trav postorder traversal
 * @param [in] start start index of from which the tree is to be constructed
 * @param [in] end end index of from which the tree is to be constructed
 * @return pointer to root of newly constructed tree
 */
std::shared_ptr<Node> constructBinaryTreeFromInPost(const std::vector<double>& in_trav, const std::vector<double>& post_trav, const int start, const int end);

std::shared_ptr<Node> constructBinaryTreeFromInPostUtil(const std::vector<double>& in_trav, const std::vector<double>& post_trav, const int start, const int end, int& post_index);
std::shared_ptr<Node> constructBinaryTreeFromPreInUtil(const std::vector<double>& pre_trav, const std::vector<double>& in_trav, const int start, const int end, int& pre_index);

/**
 * @brief Creates a binary tree of given depth and a maximum and minumum value
 * @param [in] max_depth maximum depth of tree
 * @param [in] min_val minimum value within tree
 * @param [in] max_val maximum vaue within tree
 * @return pointer to root of newly constructed tree
 */
std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, const double& min_val, const double& max_val);

/**
 * @brief Creates a full binary tree of given depth and a maximum and minumum value
 * @param [in] max_depth maximum depth of tree
 * @param [in] min_val minimum value within tree
 * @param [in] max_val maximum vaue within tree
 * @return pointer to root of newly constructed tree
 */
std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, const double& min_val, const double& max_val);

template <typename T>
std::shared_ptr<Node> generateRandomBinaryTreeUtil(const int max_depth, T& dist, std::mt19937& gen, double& random_null_thresh);

/**
 * @brief Creates a binary tree of given depth and the distribution from which the tree is to be generated
 * @param [in] max_depth maximum depth of tree
 * @param [in] dist distribution
 * @return pointer to root of newly constructed tree
 */
template<typename T>
std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, T& dist);

template <typename T>
std::shared_ptr<Node> generateRandomFullBinaryTreeUtil(const int max_depth, T& dist, std::mt19937& gen, double& random_null_thresh);

/**
 * @brief Creates a full binary tree of given depth and the distribution from which the tree is to be generated
 * @param [in] max_depth maximum depth of tree
 * @param [in] dist distribution
 * @return pointer to root of newly constructed tree
 */
template<typename T>
std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, T& dist);

// Templated function definitions:

template<typename T>
void getInOrderUtil(const std::shared_ptr<Node>& node, T& in_order) {
  if(node==nullptr) return;
  getInOrderUtil(node->left_, in_order);
  in_order.insert(in_order.end(), node->value_);
  getInOrderUtil(node->right_, in_order);
}

template<typename T>
T getInOrder(const std::shared_ptr<Node>& node) {
  T in_order;
  getInOrderUtil(node, in_order);
  return in_order;
}

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
