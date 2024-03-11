#include <gtest/gtest.h>
#include <random>
#include "geomlib/surface/polygon.hpp"

using geomlib::Polygon;

class PolygonTests: public ::testing::Test {
protected:
    std::unique_ptr<Polygon> tri;
    std::unique_ptr<Polygon> rect;
    std::unique_ptr<Polygon> penta;
    virtual void SetUp() override {
        tri = std::make_unique<Polygon>(std::vector({
            lalib::VecD<3>({ 0.0, 0.0, 0.0 }),
            lalib::VecD<3>({ 0.0, 2.0, 0.0 }),
            lalib::VecD<3>({ 1.0, 0.0, 0.0 })
        }));
        rect = std::make_unique<Polygon>(std::vector({
            lalib::VecD<3>({ 0.0, 0.0, 0.0 }),
            lalib::VecD<3>({ 0.0, 1.0, 0.0 }),
            lalib::VecD<3>({ 0.0, 1.0, 1.0 }),
            lalib::VecD<3>({ 0.0, 0.0, 1.0 })
        }));
        penta = std::make_unique<Polygon>(std::vector({
            lalib::VecD<3>({ 0.0, 0.0, 0.0 }),
            lalib::VecD<3>({ 1.0, 0.0, 0.0 }),
            lalib::VecD<3>({ 3.0, 0.5, 0.0 }),
            lalib::VecD<3>({ 1.0, 1.0, 0.0 }),
            lalib::VecD<3>({ 0.0, 1.0, 0.0 })
        }));
    }
};

TEST_F(PolygonTests, PointTest) {
    auto rand = std::mt19937(std::random_device()());
    auto rng = std::uniform_real_distribution<double>(0.0, 1.0);

    auto n_checks = 10u;
    for (auto i = 0u; i < n_checks; ++i) {
        auto s = lalib::VecD<2>({ rng(rand), rng(rand) });

        auto p_tri = tri->base() * s + tri->origin();
        EXPECT_DOUBLE_EQ(p_tri[0], tri->point(s)[0]);
        EXPECT_DOUBLE_EQ(p_tri[1], tri->point(s)[1]);
        EXPECT_DOUBLE_EQ(p_tri[2], tri->point(s)[2]);

        auto p_rect = rect->base() * s + rect->origin();
        EXPECT_DOUBLE_EQ(p_rect[0], rect->point(s)[0]);
        EXPECT_DOUBLE_EQ(p_rect[1], rect->point(s)[1]);
        EXPECT_DOUBLE_EQ(p_rect[2], rect->point(s)[2]);

        auto p_penta = penta->base() * s + penta->origin();
        EXPECT_DOUBLE_EQ(p_penta[0], penta->point(s)[0]);
        EXPECT_DOUBLE_EQ(p_penta[1], penta->point(s)[1]);
        EXPECT_DOUBLE_EQ(p_penta[2], penta->point(s)[2]);
    }
}

TEST_F(PolygonTests, NormalTest) {
    EXPECT_DOUBLE_EQ(0.0, tri->normal()[0]);
    EXPECT_DOUBLE_EQ(0.0, tri->normal()[1]);
    EXPECT_DOUBLE_EQ(-1.0, tri->normal()[2]);

    EXPECT_DOUBLE_EQ(1.0, rect->normal()[0]);
    EXPECT_DOUBLE_EQ(0.0, rect->normal()[1]);
    EXPECT_DOUBLE_EQ(0.0, rect->normal()[2]);

    EXPECT_DOUBLE_EQ(0.0, penta->normal()[0]);
    EXPECT_DOUBLE_EQ(0.0, penta->normal()[1]);
    EXPECT_DOUBLE_EQ(1.0, penta->normal()[2]);
}

TEST_F(PolygonTests, AreaTest) {
    EXPECT_DOUBLE_EQ(1.0, tri->area());
    EXPECT_DOUBLE_EQ(1.0, rect->area());
    EXPECT_DOUBLE_EQ(2.0, penta->area());
}