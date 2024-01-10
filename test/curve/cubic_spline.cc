#include "../../include/curve/cubic_spline.hpp"
#include <numbers>
#include <gtest/gtest.h>

class SplineTests: public ::testing::Test {
protected:
	std::unique_ptr<geomlib::CubicSpline<2>> pline;
	std::vector<lalib::VecD<2>> plist;

	virtual void SetUp() {
		plist = {
			circle(0.0), circle(0.25), circle(0.5), circle(0.75), circle(1.0)
		};
		pline = std::make_unique<geomlib::CubicSpline<2>>(plist);
	}

	static auto circle(double s) noexcept -> lalib::VecD<2> {
		return lalib::VecD<2> ({
			std::cos(0.5 * std::numbers::pi * s),
			std::sin(0.5 * std::numbers::pi * s),
		});
	}
};

TEST_F(SplineTests, LengthTest) {
	auto curve = geomlib::CubicSpline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	
	ASSERT_NEAR(std::sqrt(2), curve.length(), 1e-10);
}

TEST_F(SplineTests, NumOfSegmentsTest) {
	ASSERT_EQ(4, pline->n_segments());
} 

TEST_F(SplineTests, InterpolationTest) {
	EXPECT_NEAR(circle(0.0)[0], pline->point(0.0)[0], 1e-10);
	EXPECT_NEAR(circle(0.0)[1], pline->point(0.0)[1], 1e-10);

	EXPECT_NEAR(circle(0.5)[0], pline->point(0.5)[0], 1e-10);
	EXPECT_NEAR(circle(0.5)[1], pline->point(0.5)[1], 1e-10);

	EXPECT_NEAR(circle(1.0)[0], pline->point(1.0)[0], 1e-10);
	EXPECT_NEAR(circle(1.0)[1], pline->point(1.0)[1], 1e-10);
}

TEST_F(SplineTests, TangentTest) {
	auto curve = geomlib::CubicSpline<2>({
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

TEST_F(SplineTests, DerivationTest) {
	auto curve = geomlib::CubicSpline<2>({
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
