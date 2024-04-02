#pragma once
#include <functional>
#include <numeric>
#include <ranges>
#include <stack>
#include <cassert>
#include <vector>
#include <utility>
#include <memory>
#include <algorithm>

namespace geomlib::nns {

template<typename T, size_t N, typename V = T, class Proj = std::identity> class KdTree;


template<typename T, size_t N, typename V = T, class Proj = std::identity>
auto make_kd_tree(std::vector<T>&& nodes, Proj&& proj = {}) -> KdTree<T, N, V, Proj>;


template<typename T, size_t N, typename V, class Proj>
class KdTree {
public:
	class Node;
	enum class Child { Left, Right };

	KdTree(std::vector<T>&& nodes, Proj&& proj = {}) noexcept: _nodes( __construct(std::move(nodes), proj) ), _proj(std::move(proj)) {}  

	auto size() const noexcept -> size_t;

	auto nn_search(const V& query) const noexcept -> std::pair<double, std::weak_ptr<Node>>;
	auto knn_search(const V& query, size_t k = 1) const noexcept -> std::vector<std::pair<double, std::weak_ptr<Node>>>;
	auto radius_search(const V& query, double radius) const noexcept -> std::vector<std::pair<double, std::weak_ptr<Node>>>;

private:
	template<bool Sized> class __CandidateContainer;
	using __SizedCandidateContainer = __CandidateContainer<true>;
	using __UnsizedCandidateContainer = __CandidateContainer<false>;

	std::vector<std::shared_ptr<Node>> _nodes;
	Proj _proj;

	auto __get_entry_point() const noexcept -> std::weak_ptr<Node> { return this->_nodes.back(); }

	static auto __construct(std::vector<T>&& nodes, const Proj& proj) noexcept -> std::vector<std::shared_ptr<Node>>;
	static auto __distance(std::weak_ptr<Node> node, const V& query, Proj proj) noexcept -> double;
	static auto __get_leaf(std::weak_ptr<Node> node, const V& query, std::stack<std::pair<Child, std::weak_ptr<Node>>>& stack) noexcept -> std::weak_ptr<Node>;
};


template<typename T, size_t N, typename V, class Proj>
class KdTree<T, N, V, Proj>::Node {
public:
	Node(size_t axis, T&& value, Proj proj) noexcept: _value(std::move(value)), _axis_no(axis), _proj(proj) {}
	Node(size_t axis, T&& value, std::weak_ptr<Node> left, std::weak_ptr<Node> right, Proj proj) noexcept:
		_value(std::move(value)), _l(left), _r(right), _axis_no(axis), _proj(proj)
	{}

	auto is_leaf() const noexcept -> bool;

	auto value() const noexcept -> const T& { return this->_value; }
	auto child(Child ch) const noexcept -> std::weak_ptr<Node> { return ch == Child::Left ? this->_l : this->_r; }
	auto left() const noexcept -> std::weak_ptr<Node> { return this->_l; }
	auto right() const noexcept -> std::weak_ptr<Node>{ return this->_r; }

	auto select_lr(const V& query) const noexcept -> std::pair<Child, std::weak_ptr<Node>>;
	auto distance(const V& query) const noexcept -> double;
	auto distance_to_plane(const V& query) const noexcept -> double;

private:
	T 					_value;
	std::weak_ptr<Node> _l;
	std::weak_ptr<Node> _r;
	size_t 				_axis_no;
	Proj				_proj;
};


template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
class KdTree<T, N, V, Proj>::__CandidateContainer {
public:
	__CandidateContainer() noexcept;
	__CandidateContainer(size_t max_cands) noexcept;

	auto size() const noexcept -> size_t;
	void add(double dist, std::weak_ptr<Node> node) noexcept;

	auto nearest_distance() const noexcept -> double;
	auto furthest_distance() const noexcept -> double;

	auto get_nearest() const noexcept -> std::vector<std::pair<double, std::weak_ptr<Node>>>; 

private:
	size_t _max_cands;
	std::vector<std::pair<double, std::weak_ptr<Node>>> _cands;
};


template<typename T, size_t N, typename V, class Proj>
auto make_kd_tree(std::vector<T>&& nodes, Proj&& proj) -> KdTree<T, N, V, Proj> {
	return KdTree<T, N, V, Proj>(std::move(nodes), std::move(proj));
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::nn_search(const V& query) const noexcept -> std::pair<double, std::weak_ptr<Node>> {
	auto root = this->__get_entry_point();
	auto stack = std::stack<std::pair<Child, std::weak_ptr<Node>>>();

	auto leaf = __get_leaf(root, query, stack);
	auto d = __distance(leaf, query, this->_proj);
	auto nearest = std::make_pair(d, leaf);

	while (!stack.empty()) {
		auto [lr, parent] = stack.top();
		stack.pop();

		{
			auto tmp_d = __distance(parent, query, this->_proj);
			if (tmp_d < d) { 
				d = tmp_d; 
				nearest = std::make_pair(tmp_d, parent);
			}
		}

		if (parent.lock()->distance_to_plane(query) < d) {
			auto sub_node = parent.lock()->child(lr);
			if (!sub_node.expired()) {
				leaf = __get_leaf(sub_node, query, stack);				
				auto tmp_d = __distance(leaf, query, this->_proj);
				if (tmp_d < d) { 
					d = tmp_d; 
					nearest = std::make_pair(tmp_d, leaf);
				} 
			}
		}

	}

	return nearest;
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::knn_search(const V& query, size_t k) const noexcept -> std::vector<std::pair<double, std::weak_ptr<Node>>> {
	auto cands = __SizedCandidateContainer(k);
	auto root = this->__get_entry_point();
	auto stack = std::stack<std::pair<Child, std::weak_ptr<Node>>>();

	auto leaf = __get_leaf(root, query, stack);
	auto d = __distance(leaf, query, this->_proj);
	cands.add(d, leaf);

	while (!stack.empty()) {
		auto [lr, parent] = stack.top();
		stack.pop();

		if (leaf.lock() != parent.lock()) {
			auto tmp_d = __distance(parent, query, this->_proj);
			cands.add(tmp_d, parent);
		}

		if (parent.lock()->distance_to_plane(query) < cands.furthest_distance()) {
			auto sub_node = parent.lock()->child(lr);
			if (!sub_node.expired()) {
				leaf = __get_leaf(sub_node, query, stack);				

				auto tmp_d = __distance(leaf, query, this->_proj);
				cands.add(tmp_d, leaf);
			}
		}
	}
	return cands.get_nearest();
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::radius_search(const V& query, double radius) const noexcept -> std::vector<std::pair<double, std::weak_ptr<Node>>> {
	auto cands = __UnsizedCandidateContainer();
	auto root = this->__get_entry_point();
	auto stack = std::stack<std::pair<Child, std::weak_ptr<Node>>>();

	auto leaf = __get_leaf(root, query, stack);
	auto d = __distance(leaf, query, this->_proj);
	if (d < radius) { cands.add(d, leaf); } 

	while (!stack.empty()) {
		auto [lr, parent] = stack.top();
		stack.pop();

		if (leaf.lock() != parent.lock()) {
			auto tmp_d = __distance(parent, query, this->_proj);
			if (tmp_d < radius) { cands.add(tmp_d, parent); }
		}

		if (parent.lock()->distance_to_plane(query) < radius) { 
			auto sub_node = parent.lock()->child(lr);
			if (!sub_node.expired()) {
				leaf = __get_leaf(sub_node, query, stack);				

				auto tmp_d = __distance(leaf, query, this->_proj);
				if (tmp_d < radius) { cands.add(tmp_d, leaf); }
			}
		}
	}
	return cands.get_nearest();
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::__construct(std::vector<T>&& values, const Proj& proj) noexcept -> std::vector<std::shared_ptr<Node>> {
	auto nodes = std::vector<std::shared_ptr<Node>>();
	auto indices = std::vector<size_t>(values.size(), 0);
	std::iota(indices.begin(), indices.end(), 0);

	// Push to stack with depth-first search
	auto index_stack = std::stack<std::pair<size_t, std::vector<size_t>>>();
	auto node_stack = std::stack<std::pair<size_t, std::weak_ptr<Node>>>();
	auto value_stack = std::stack<size_t>();

	size_t depth = 0u;
	index_stack.emplace( std::make_pair(depth, std::move(indices)) );
	while (!index_stack.empty() or !node_stack.empty()) {
		if ( !index_stack.empty() and (node_stack.empty() or index_stack.top().first >= node_stack.top().first)) {
			auto cur_indices = index_stack.top().second;
			index_stack.pop();

			// Create leaf node
			if (cur_indices.size() == 0) {
				auto cur_depth = depth;
				if (!node_stack.empty() and node_stack.top().first == depth) { --depth; }
				node_stack.emplace( std::make_pair(cur_depth, std::weak_ptr<Node>()) );
			}
			else if (cur_indices.size() == 1) {
				auto cur_depth = depth;
				auto new_node = nodes.emplace_back( 
					std::make_shared<Node>( cur_depth % N, std::move(values[cur_indices[0]]), proj )
				);
				if (!node_stack.empty() and node_stack.top().first == depth) { --depth; }
				node_stack.emplace( std::make_pair( cur_depth, new_node));
			}
			// Push the next-depth state to the stacks
			else {
				size_t mid = cur_indices.size() / 2;
				auto mid_iter = cur_indices.begin() + mid;
				std::ranges::nth_element(cur_indices, mid_iter, [&](const auto& a, const auto& b) { 
					auto axis = depth % N;
					return std::invoke(proj, values[a])[axis] < std::invoke(proj, values[b])[axis];
				});

				++depth;
				index_stack.emplace( std::make_pair(depth, std::vector(cur_indices.begin(), cur_indices.begin() + mid)) );
				index_stack.emplace( std::make_pair(depth, std::vector(cur_indices.begin() + mid + 1, cur_indices.end())) );
				value_stack.emplace( cur_indices[mid] );
			}
		}
		// Create internal node
		else {
			auto l = node_stack.top().second;
			node_stack.pop();
			auto r = node_stack.top().second;
			node_stack.pop();

			auto new_node = nodes.emplace_back( 
				std::make_shared<Node>( depth % N, std::move(values[value_stack.top()]), l, r, proj )
			);
			value_stack.pop();

			auto cur_depth = depth;
			if (!node_stack.empty() and node_stack.top().first == depth) { --depth; }
			if (cur_depth != 0) {
				node_stack.emplace( std::make_pair( cur_depth, new_node));
			}
		}
	}
	return nodes;
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::__distance(std::weak_ptr<Node> node, const V& query, Proj proj) noexcept -> double {
	return ( std::invoke(proj, node.lock()->value()) - query ).norm2();
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::__get_leaf(std::weak_ptr<Node> node, const V& query, std::stack<std::pair<Child, std::weak_ptr<Node>>>& stack) noexcept -> std::weak_ptr<Node> {
	// Depth-first search
	while (!node.lock()->is_leaf()) { 
		auto [lr, new_node] = node.lock()->select_lr(query);

		stack.push(std::make_pair(
			lr == Child::Left ? Child::Right: Child::Left, 
			node
		));

		if (new_node.expired()) { return node; }
		else { node = new_node; }
	}
	return node;	
}



// ==== Implementation of KdTree::Node ===== //

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::Node::is_leaf() const noexcept -> bool { 
	return this->_l.expired() and this->_r.expired();
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::Node::select_lr(const V& query) const noexcept -> std::pair<Child, std::weak_ptr<Node>> { 
	auto val = std::invoke(this->_proj, this->_value)[this->_axis_no];
	if ( query[this->_axis_no] < val ) {
		return std::make_pair(Child::Left, this->_l);
	} else {
		return std::make_pair(Child::Right, this->_r);
	}
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::Node::distance(const V& query) const noexcept -> double { 
	return ( std::invoke(this->_proj, this->_value) - query ).norm();
}

template<typename T, size_t N, typename V, class Proj>
auto KdTree<T, N, V, Proj>::Node::distance_to_plane(const V& query) const noexcept -> double { 
	return std::abs(	
		std::invoke(this->_proj, this->_value)[this->_axis_no] - query[this->_axis_no]
	);
}


// ==== Implementation of KdTree::__CandidateContainer ===== //

template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
KdTree<T, N, V, Proj>::__CandidateContainer<Sized>::__CandidateContainer() noexcept: 
	_max_cands(std::numeric_limits<size_t>::max()),
	_cands()
{  
	if constexpr (Sized) { static_assert([]{ return false; }()); }
}

template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
KdTree<T, N, V, Proj>::__CandidateContainer<Sized>::__CandidateContainer(size_t max_cands) noexcept: 
	_max_cands(Sized ? max_cands : std::numeric_limits<size_t>::max()),
	_cands()
{  
	_cands.reserve(max_cands);
}

template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
auto KdTree<T, N, V, Proj>::__CandidateContainer<Sized>::size() const noexcept -> size_t { 
	return this->_cands.size();
}

template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
void KdTree<T, N, V, Proj>::__CandidateContainer<Sized>::add(double dist, std::weak_ptr<Node> node) noexcept { 
	if constexpr (Sized) {
		if (this->_cands.size() == this->_max_cands) {
			if (dist < this->_cands.back().first) 	{ this->_cands.pop_back(); }
			else 									{ return; }
		}
	}
	auto insert_pos = std::ranges::upper_bound(this->_cands, dist, {}, &std::pair<double, std::weak_ptr<Node>>::first);
	this->_cands.insert(insert_pos, std::make_pair(dist, node));
}

template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
auto KdTree<T, N, V, Proj>::__CandidateContainer<Sized>::nearest_distance() const noexcept -> double { 
	assert(!this->_cands.empty());
	return this->_cands.front().first;
}
template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
auto KdTree<T, N, V, Proj>::__CandidateContainer<Sized>::furthest_distance() const noexcept -> double { 
	auto idx = this->_cands.size() - 1;
	return this->_cands[idx].first;
}

template<typename T, size_t N, typename V, class Proj>
template<bool Sized>
auto KdTree<T, N, V, Proj>::__CandidateContainer<Sized>::get_nearest() const noexcept -> std::vector<std::pair<double, std::weak_ptr<Node>>> {  
	auto cands_view = std::ranges::views::all(this->_cands);
	return std::vector(cands_view.begin(), cands_view.end());
}

}
