#include <gtest/gtest.h>
#include <functional>
#include "../include/curve.hpp"
#include "../include/affine.hpp"

TEST(CurveCommonTests, ApproxLengthTest) {
	auto seg = geomlib::Segment<2>(lalib::VecD<2>({ 0.0, 0.0 }), lalib::VecD<2>({ 1.0, 1.0 }));
	ASSERT_NEAR(std::sqrt(2), geomlib::approx_length(0.0, 1.0, seg, 10), 1e-10);

	auto polyline = geomlib::CubicSpline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	ASSERT_NEAR(std::sqrt(2), geomlib::approx_length(0.0, 1.0, polyline, 100), 1e-6); 

	auto spline = geomlib::CubicSpline<2>({
		lalib::VecD<2>({ 0.0, 0.0 }),
		lalib::VecD<2>({ 0.2, 0.2 }),
		lalib::VecD<2>({ 0.7, 0.7 }),
		lalib::VecD<2>({ 1.0, 1.0 })
	});
	ASSERT_NEAR(std::sqrt(2), geomlib::approx_length(0.0, 1.0, spline, 100), 1e-6); 
}

TEST(CurveCommonTests, ResamplingTest) {
    auto plist = std::vector {
        lalib::VecD<2>({ 0.0, 0.0 }),
        lalib::VecD<2>({ 0.2, 0.0 }),
        lalib::VecD<2>({ 0.5, 0.0 }),
        lalib::VecD<2>({ 1.0, 0.0 }),
    };
    auto spline = geomlib::CubicSpline<2>(plist);
    auto sampled = geomlib::sample(spline, 11, std::identity(), 1e-6);

    for (auto&& seg: sampled) {
        EXPECT_NEAR(seg.length(), 0.1, 1e-6);
    }
}

TEST(CurveCommonTests, ApproxResamplingTest) {
    auto plist = std::vector {
        lalib::VecD<2>({ 0.0, 0.0 }),
        lalib::VecD<2>({ 0.2, 0.0 }),
        lalib::VecD<2>({ 0.5, 0.0 }),
        lalib::VecD<2>({ 1.0, 0.0 }),
    };
    auto spline = geomlib::CubicSpline<2>(plist);
    auto affine = geomlib::rotate2d(0.0, lalib::VecD<2>::filled(0.0));
    auto curve = geomlib::transform_lazy(std::move(spline), std::move(affine));

    auto sampled = geomlib::sample(curve, 11, std::identity(), 1e-6);

    for (auto&& seg: sampled) {
        EXPECT_NEAR(seg.length(), 0.1, 1e-6);
    }
}

TEST(CurveCommonTests, NearestPointTest) {
    auto plist = std::vector {
        lalib::VecD<2>({ 1.0, 0.0 }),
        lalib::VecD<2>({ 0.8, -0.1 }),
        lalib::VecD<2>({ 0.2, -0.1 }),
        lalib::VecD<2>({ 0.0, 0.0 }),
        lalib::VecD<2>({ 0.2, 0.1 }),
        lalib::VecD<2>({ 0.8, 0.1 }),
        lalib::VecD<2>({ 1.0, 0.01 })
    };
    auto curve = geomlib::CubicSpline(std::move(plist));

    {
        auto exact_s = 0.7;
        auto query = curve.point(exact_s);
        auto found_s = geomlib::nearest_on(curve, query, 10);
        EXPECT_NEAR(exact_s, found_s, 1e-6);
    }

    {
        auto exact_s = 0.0;
        auto query = curve.point(exact_s);
        auto found_s = geomlib::nearest_on(curve, query, 10);
        EXPECT_NEAR(exact_s, found_s, 1e-6);
    }

    {
        auto exact_s = 1.0;
        auto query = curve.point(exact_s);
        auto found_s = geomlib::nearest_on(curve, query, 10);
        EXPECT_NEAR(exact_s, found_s, 1e-6);
    }

    {
        // Off the curve
        auto query_edge2 = lalib::VecD<2>({ 2.0, 0.01 });
        auto found_s_edge2 = geomlib::nearest_on(curve, query_edge2, 10);
        EXPECT_NEAR(1.0, found_s_edge2, 1e-6);
    }
}