#include <gtest/gtest.h>
#include <ranges>
#include "../../include/curve/polyline.hpp"

static_assert(std::input_iterator<geomlib::Polyline<2>::SegmentViewIterator>);
static_assert(std::sentinel_for<geomlib::Polyline<2>::SegmentViewIterator, geomlib::Polyline<2>::SegmentViewIterator>);
static_assert(std::ranges::input_range<geomlib::Polyline<2>>);

TEST(PolylineTests, LengthTest) {
	auto curve = geomlib::Polyline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	
	ASSERT_NEAR(std::sqrt(2), curve.length(), 1e-10);
}

TEST(PolylineTests, NumOfSegmentTest) {
	auto curve = geomlib::Polyline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	
	ASSERT_EQ(3, curve.n_segments());
}

TEST(PolylineTests, SegmentIterationTest) {
	auto plist = std::vector({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	auto curve = geomlib::Polyline<2>(plist);

	auto i = 0u;
	for (auto&& seg: curve) {
		ASSERT_DOUBLE_EQ(seg.point(0.0)[0], plist[i][0]);
		ASSERT_DOUBLE_EQ(seg.point(0.0)[1], plist[i][1]);

		ASSERT_DOUBLE_EQ(seg.point(1.0)[0], plist[i + 1][0]);
		ASSERT_DOUBLE_EQ(seg.point(1.0)[1], plist[i + 1][1]);

		++i;
	}
}

TEST(PolylineTests, InterPolationTest) {
	auto curve = geomlib::Polyline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	
    ASSERT_DOUBLE_EQ(0.0, curve.point(0.0)[0]);
    ASSERT_DOUBLE_EQ(0.0, curve.point(0.0)[1]);

    ASSERT_DOUBLE_EQ(0.5, curve.point(0.5)[0]);
    ASSERT_DOUBLE_EQ(0.5, curve.point(0.5)[1]);

    ASSERT_DOUBLE_EQ(1.0, curve.point(1.0)[0]);
    ASSERT_DOUBLE_EQ(1.0, curve.point(1.0)[1]);
}

TEST(PolylineTests, TangentTest) {
	auto curve = geomlib::Polyline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});

	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.0)[0]);
	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.0)[1]);

	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.1)[0]);
	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.1)[1]);

	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.5)[0]);
	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.5)[1]);

	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.7)[0]);
	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.7)[1]);

	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.9)[0]);
	EXPECT_DOUBLE_EQ(1.0 / std::sqrt(2), curve.tangent(0.9)[1]);
}

TEST(PolylineTests, DerivationTest) {
	auto curve = geomlib::Polyline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});

	EXPECT_DOUBLE_EQ(0.2 / sqrt(0.2*0.2 + 0.2*0.2), curve.tangent(0.0)[0]);
	EXPECT_DOUBLE_EQ(0.2 / sqrt(0.2*0.2 + 0.2*0.2), curve.tangent(0.0)[1]);

	EXPECT_DOUBLE_EQ(0.5 / (sqrt(0.7*0.7 + 0.7*0.7) - sqrt(0.2*0.2 + 0.2*0.2)), curve.tangent(0.5)[0]);
	EXPECT_DOUBLE_EQ(0.5 / (sqrt(0.7*0.7 + 0.7*0.7) - sqrt(0.2*0.2 + 0.2*0.2)), curve.tangent(0.5)[1]);
}
