#include <memory>
#include <vector>

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
std::shared_ptr<Node> generateRandomBinaryTree(const int max_depth, const double min_val, const double max_val);
std::shared_ptr<Node> generateRandomFullBinaryTree(const int max_depth, const double min_val, const double max_val);
}; // namespace binary_tree
