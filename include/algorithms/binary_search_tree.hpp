#ifndef BINARY_SEARCH_TREE_HPP
#define BINARY_SEARCH_TREE_HPP

#include "algorithms/binary_tree.hpp"

namespace binary_search_tree {

std::shared_ptr<binary_tree::Node> constructTreeFromPreOrderUtil_N2(const std::vector<int>& pre_order, int start, int end, int& pre_index);
std::shared_ptr<binary_tree::Node> constructTreeFromPreOrderUtil_N(const std::vector<double>& pre_order, double min, double max, int& pre_index);
std::shared_ptr<binary_tree::Node> constructTreeFromPreOrderUtilStack_N(const std::vector<double>& pre_order);
std::shared_ptr<binary_tree::Node> constructTreeFromPreOrder(const std::vector<double>& pre_order);

std::shared_ptr<binary_tree::Node> insertIntoBST(std::shared_ptr<binary_tree::Node> node, const double value);
std::tuple<bool, std::shared_ptr<binary_tree::Node>> findInBST(std::shared_ptr<binary_tree::Node> node, const double value);
bool deleteFromBST(std::shared_ptr<binary_tree::Node>& node, const double value);

std::shared_ptr<binary_tree::Node> findLowestCommonAncestorUtil(const std::shared_ptr<binary_tree::Node>& node, double v1, double v2);
std::shared_ptr<binary_tree::Node> findLowestCommonAncestor(const std::shared_ptr<binary_tree::Node>& node, double v1, double v2);

std::tuple<bool, std::vector<std::pair<double, double>>> checkPairSumBST(const std::shared_ptr<binary_tree::Node> node, const double value);
}; // namespace binary_search_tree

#endif  // BINARY_SEARCH_TREE_HPP
