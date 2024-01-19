#include <gtest/gtest.h>
#include <functional>
#include "../include/curve.hpp"

TEST(CurveCommonTests, ResamplingTest) {
    auto plist = std::vector {
        lalib::VecD<2>({ 0.0, 0.0 }),
        lalib::VecD<2>({ 0.2, 0.0 }),
        lalib::VecD<2>({ 0.5, 0.0 }),
        lalib::VecD<2>({ 1.0, 0.0 }),
    };
    auto spline = geomlib::CubicSpline<2>(plist);
    auto sampled = geomlib::sample(spline, 11, std::identity());

    for (auto&& seg: sampled) {
        EXPECT_NEAR(seg.length(), 0.1, 1e-6);
    }
}