#include <gtest/gtest.h>
#include "../../include/curve/polyline.hpp"

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
