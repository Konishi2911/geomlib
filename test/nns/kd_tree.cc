#include "geomlib/nns/kd_tree.hpp"
#include "lalib//vec.hpp"
#include <gtest/gtest.h>
#include <random>
#include <iostream>

class KdTreeTests: public ::testing::Test {
protected:
	using V3 = lalib::VecD<3>;

	static constexpr size_t n = 100;
	std::unique_ptr<geomlib::nns::KdTree<V3, 3>> 	kd_tree;
	std::vector<V3> 								points;
	std::vector<std::pair<size_t, V3>> 				pair_points;

	KdTreeTests(): kd_tree(), points() {
		auto d = std::random_device()();
		std::cout << " Seed = " << d << std::endl;
		//auto rng = std::mt19937(483615696);
		auto rng = std::mt19937(d);
		auto dist = std::uniform_real_distribution<double>(-10.0, 10.0);

		for (auto i = 0u; i < n; ++i) {
			auto&p = points.emplace_back(
				V3({ dist(rng), dist(rng), dist(rng) })
			);
			pair_points.emplace_back( std::make_pair( i, p ));
		}
		
		auto tmp = points;
		kd_tree = std::make_unique<geomlib::nns::KdTree<V3, 3, V3>>(std::move(tmp));

		auto tree = geomlib::nns::make_kd_tree<std::pair<size_t, V3>, 3, V3>(std::vector(pair_points), &std::pair<size_t, V3>::second);
	}

	virtual void SetUp() override {
	}

};

TEST_F(KdTreeTests, NNTest) {
	auto query = lalib::SizedVec<double, 3>({ 0.0, 0.0, 0.0 });
	auto tmp = points;
	std::ranges::sort(tmp.begin(), tmp.end(), [&](auto&& a, auto&& b) { return (query - a).norm2() < (query - b).norm2(); });
	
	auto nn_result = kd_tree->nn_search(query);

	ASSERT_DOUBLE_EQ((nn_result.second.lock()->value() - query).norm2(), (tmp[0] - query).norm2());
	ASSERT_DOUBLE_EQ(nn_result.second.lock()->value()[0], tmp[0][0]);
	ASSERT_DOUBLE_EQ(nn_result.second.lock()->value()[1], tmp[0][1]);
	ASSERT_DOUBLE_EQ(nn_result.second.lock()->value()[2], tmp[0][2]);
}

TEST_F(KdTreeTests, kNNTest) {
	auto query = lalib::SizedVec<double, 3>({ 0.0, 0.0, 0.0 });
	auto tmp = points;
	std::ranges::sort(tmp.begin(), tmp.end(), [&](auto&& a, auto&& b) { return (query - a).norm2() < (query - b).norm2(); });
	
	auto knn_result = kd_tree->knn_search(query, 5);

	ASSERT_DOUBLE_EQ((knn_result[0].second.lock()->value() - query).norm2(), (tmp[0] - query).norm2());
	ASSERT_DOUBLE_EQ(knn_result[0].second.lock()->value()[0], tmp[0][0]);
	ASSERT_DOUBLE_EQ(knn_result[0].second.lock()->value()[1], tmp[0][1]);
	ASSERT_DOUBLE_EQ(knn_result[0].second.lock()->value()[2], tmp[0][2]);

	ASSERT_DOUBLE_EQ((knn_result[1].second.lock()->value() - query).norm2(), (tmp[1] - query).norm2());
	ASSERT_DOUBLE_EQ((knn_result[2].second.lock()->value() - query).norm2(), (tmp[2] - query).norm2());
}
