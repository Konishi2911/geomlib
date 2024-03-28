#include <gtest/gtest.h>
#include "geomlib/curve/polynomial.hpp"
#include "geomlib/curve/curve_concepts.hpp"

static_assert(geomlib::Curve<geomlib::PolynomialCurve<3>>);

TEST(PolynomialCurveTests, PointTest) {
    // quadratic function
    auto curve = geomlib::PolynomialCurve<3>(std::vector {
        lalib::VecD<3>({ 0.0, 0.0, 0.0 }),
        lalib::VecD<3>({ 1.0, 1.0, 0.0 }),
        lalib::VecD<3>({ 0.0, 0.0, 1.0 })
    });

    ASSERT_DOUBLE_EQ(0.0, curve.point(0.0)[0]);
    ASSERT_DOUBLE_EQ(0.0, curve.point(0.0)[1]);
    ASSERT_DOUBLE_EQ(0.0, curve.point(0.0)[2]);

    ASSERT_DOUBLE_EQ(0.5, curve.point(0.5)[0]);
    ASSERT_DOUBLE_EQ(0.5, curve.point(0.5)[1]);
    ASSERT_DOUBLE_EQ(0.25, curve.point(0.5)[2]);

    ASSERT_DOUBLE_EQ(1.0, curve.point(1.0)[0]);
    ASSERT_DOUBLE_EQ(1.0, curve.point(1.0)[1]);
    ASSERT_DOUBLE_EQ(1.0, curve.point(1.0)[2]);
}

TEST(PolynomialCurveTests, DerivationTest) {
    // Pure quadratic function
    auto curve = geomlib::PolynomialCurve<3>(std::vector {
        lalib::VecD<3>({ 0.0, 0.0, 0.0 }),
        lalib::VecD<3>({ 1.0, 1.0, 0.0 }),
        lalib::VecD<3>({ 0.0, 0.0, 1.0 })
    });

    ASSERT_DOUBLE_EQ(1.0, curve.deriv(0.0)[0]);
    ASSERT_DOUBLE_EQ(1.0, curve.deriv(0.0)[1]);
    ASSERT_DOUBLE_EQ(0.0, curve.deriv(0.0)[2]);

    ASSERT_DOUBLE_EQ(1.0, curve.deriv(0.5)[0]);
    ASSERT_DOUBLE_EQ(1.0, curve.deriv(0.5)[1]);
    ASSERT_DOUBLE_EQ(1.0, curve.deriv(0.5)[2]);

    ASSERT_DOUBLE_EQ(1.0, curve.deriv(1.0)[0]);
    ASSERT_DOUBLE_EQ(1.0, curve.deriv(1.0)[1]);
    ASSERT_DOUBLE_EQ(2.0, curve.deriv(1.0)[2]);
}