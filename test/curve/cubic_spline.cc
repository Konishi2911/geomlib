#include "../../include/curve/cubic_spline.hpp"
#include "../../include/affine.hpp"
#include <iostream>
#include <numbers>
#include <gtest/gtest.h>

template<size_t N>
void print_spline(const geomlib::CubicSpline<N>& spline, size_t n) noexcept {
	for (auto i = 0u; i <= n; ++i) {
		auto denom = static_cast<double>(n);
		std::cout << spline.point(i / denom)[0] << ", " << spline.point(i / denom)[1] << std::endl;
	}
}

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

TEST_F(SplineTests, InterLengthTest) {
	auto curve = geomlib::CubicSpline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.0 }),
		lalib::VecD<2>({ 0.4, 0.0 }),
		lalib::VecD<2>({ 0.6, 0.0 }),
		lalib::VecD<2>({ 0.8, 0.0 }),
		lalib::VecD<2>({ 1.0, 0.0 })
	});
	
	ASSERT_NEAR(0.1, curve.length(0.1), 1e-10);
	ASSERT_NEAR(0.7, curve.length(0.7), 1e-10);
}

TEST_F(SplineTests, NumOfSegmentsTest) {
	ASSERT_EQ(4, pline->n_segments());
} 

TEST_F(SplineTests, InterpolationTest) {
	// Equal-spacing circle
	EXPECT_NEAR(circle(0.0)[0], pline->point(0.0)[0], 1e-10);
	EXPECT_NEAR(circle(0.0)[1], pline->point(0.0)[1], 1e-10);

	EXPECT_NEAR(circle(0.5)[0], pline->point(0.5)[0], 1e-10);
	EXPECT_NEAR(circle(0.5)[1], pline->point(0.5)[1], 1e-10);

	EXPECT_NEAR(circle(1.0)[0], pline->point(1.0)[0], 1e-10);
	EXPECT_NEAR(circle(1.0)[1], pline->point(1.0)[1], 1e-10);

	// Unequal-spacing line
	auto plist = std::vector<lalib::VecD<2>>({
        lalib::VecD<2>({ 0.0, 0.0 }),
        lalib::VecD<2>({ 0.5, 0.0 }),
        lalib::VecD<2>({ 0.7, 0.0 }),
        lalib::VecD<2>({ 1.0, 0.0 })
    });
	auto spline = geomlib::CubicSpline<2>(plist);

	//print_spline(spline, 50);
	
    ASSERT_DOUBLE_EQ(0.0, spline.point(0.0)[0]);
    ASSERT_DOUBLE_EQ(0.0, spline.point(0.0)[1]);

    ASSERT_DOUBLE_EQ(1.0, spline.point(1.0)[0]);
    ASSERT_DOUBLE_EQ(0.0, spline.point(1.0)[1]);	
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

TEST_F(SplineTests, AffineTransformationTest) {
	auto curve = geomlib::CubicSpline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	
	auto affine = geomlib::rotate2d(std::numbers::pi / 4.0, lalib::VecD<2>({ 1.0, 1.0 }));
	auto half_scale = geomlib::Affine(lalib::MatD<2, 2>({
		0.5, 0.0,
		0.0, 0.5
	}), lalib::VecD<2>::filled(0.0));

	static_assert(geomlib::AffineTransformableCurve<decltype(curve), 2>);

	curve.transform(affine);
	ASSERT_NEAR(std::sqrt(2), curve.length(), 1e-10);

	auto half_curve = curve.transformed(half_scale);
	ASSERT_NEAR(std::sqrt(2) / 2.0, half_curve.length(), 1e-10);
}