#include <functional>
#include <gtest/gtest.h>
#include <iterator>
#include <random>
#include "geomlib/surface/polygon.hpp"

using geomlib::Polygon;
using geomlib::PolygonView;

class PolygonTests: public ::testing::Test {
protected:
    std::unique_ptr<Polygon> tri;
    std::unique_ptr<Polygon> rect;
    std::unique_ptr<Polygon> rect2;
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
        rect2 = std::make_unique<Polygon>(std::vector({
            lalib::VecD<3>({ 0.0, 1.0, 0.0 }),
            lalib::VecD<3>({ 0.0, 1.0, 1.0 }),
            lalib::VecD<3>({ 0.0, 0.0, 1.0 }),
            lalib::VecD<3>({ 0.0, 0.0, 0.0 })
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

class PolygonViewTests: public ::testing::Test {
protected:
    std::vector<lalib::VecD<3>> points = std::vector({
        lalib::VecD<3>({0.0, 0.0, 0.0}),
        lalib::VecD<3>({0.0, 2.0, 0.0}),
        lalib::VecD<3>({1.0, 0.0, 0.0}),
        lalib::VecD<3>({0.0, 1.0, 0.0}),
        lalib::VecD<3>({0.0, 1.0, 1.0}),
        lalib::VecD<3>({0.0, 0.0, 1.0}),
        lalib::VecD<3>({3.0, 0.5, 0.0}),
        lalib::VecD<3>({1.0, 0.0, 0.0}),
        lalib::VecD<3>({1.0, 1.0, 0.0})
    });
    std::unique_ptr<PolygonView> tri;
    std::unique_ptr<PolygonView> rect;
    std::unique_ptr<PolygonView> rect2;
    std::unique_ptr<PolygonView> penta;
    virtual void SetUp() override {
        tri = std::make_unique<PolygonView>(std::vector({
            std::cref(points[0]),
            std::cref(points[1]),
            std::cref(points[2])
        }));
        rect = std::make_unique<PolygonView>(std::vector({
            std::cref(points[0]),
            std::cref(points[3]),
            std::cref(points[4]),
            std::cref(points[5])
        }));
        rect2 = std::make_unique<PolygonView>(std::vector({
            std::cref(points[3]),
            std::cref(points[4]),
            std::cref(points[5]),
            std::cref(points[0])
        }));
        penta = std::make_unique<PolygonView>(std::vector({
            std::cref(points[0]),
            std::cref(points[2]),
            std::cref(points[6]),
            std::cref(points[8]),
            std::cref(points[3])
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

TEST_F(PolygonTests, FootTest) {
    auto query1 = lalib::VecD<3>({ 1.0, 0.5, 0.5 });
    auto foot1 = rect->foot(query1);

    auto query2 = lalib::VecD<3>({ 0.2, 0.8, 2.0 });
    auto foot2 = tri->foot(query2);

    EXPECT_DOUBLE_EQ(0.0, foot1[0]);
    EXPECT_DOUBLE_EQ(0.5, foot1[1]);
    EXPECT_DOUBLE_EQ(0.5, foot1[2]);

    EXPECT_DOUBLE_EQ(0.2, foot2[0]);
    EXPECT_DOUBLE_EQ(0.8, foot2[1]);
    EXPECT_DOUBLE_EQ(0.0, foot2[2]);
}

TEST_F(PolygonTests, IntersectionTest) {
    auto query = lalib::VecD<3>({ 0.0, 0.0, 1.0 });
    auto v = lalib::VecD<3>({ 0.0, 1.0 / std::sqrt(2), -1.0 / std::sqrt(2) });
    auto intersec = tri->intersection(query, v);

    EXPECT_DOUBLE_EQ(0.0, intersec[0]);
    EXPECT_DOUBLE_EQ(1.0, intersec[1]);
    EXPECT_DOUBLE_EQ(0.0, intersec[2]);
}

TEST_F(PolygonTests, LocalConversionTest) {
    auto query1 = lalib::VecD<3>({ 0.0, 0.5, 0.5 });
    auto local1 = rect->local(query1);
    auto restore1 = rect->point(local1);

    EXPECT_DOUBLE_EQ(restore1[0], query1[0]);
    EXPECT_DOUBLE_EQ(restore1[1], query1[1]);
    EXPECT_DOUBLE_EQ(restore1[2], query1[2]);

    auto query2 = lalib::VecD<3>({ 0.0, 0.2, 0.2 });
    auto local2 = rect2->local(query2);
    auto restore2 = rect2->point(local2);

    EXPECT_DOUBLE_EQ(restore2[0], query2[0]);
    EXPECT_DOUBLE_EQ(restore2[1], query2[1]);
    EXPECT_DOUBLE_EQ(restore2[2], query2[2]);
}

TEST_F(PolygonTests, InclusionTest) {
    auto query1 = lalib::VecD<3>({ 0.0, 0.5, 0.5 });
    auto query2 = lalib::VecD<3>({ 0.0, 2.0, 0.5 });

    EXPECT_TRUE(rect->check_inclusion(query1));
    EXPECT_FALSE(rect->check_inclusion(query2));
}

TEST_F(PolygonViewTests, PointTest) {
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

TEST_F(PolygonViewTests, NormalTest) {
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

TEST_F(PolygonViewTests, AreaTest) {
    EXPECT_DOUBLE_EQ(1.0, tri->area());
    EXPECT_DOUBLE_EQ(1.0, rect->area());
    EXPECT_DOUBLE_EQ(2.0, penta->area());
}

TEST_F(PolygonViewTests, FootTest) {
    auto query = lalib::VecD<3>({ 1.0, 0.5, 0.5 });
    auto foot = rect->foot(query);

    EXPECT_DOUBLE_EQ(0.0, foot[0]);
    EXPECT_DOUBLE_EQ(0.5, foot[1]);
    EXPECT_DOUBLE_EQ(0.5, foot[2]);
}

TEST_F(PolygonViewTests, IntersectionTest) {
    auto query = lalib::VecD<3>({ 0.0, 0.0, 1.0 });
    auto v = lalib::VecD<3>({ 0.0, 1.0 / std::sqrt(2), -1.0 / std::sqrt(2) });
    auto intersec = tri->intersection(query, v);

    EXPECT_DOUBLE_EQ(0.0, intersec[0]);
    EXPECT_DOUBLE_EQ(1.0, intersec[1]);
    EXPECT_DOUBLE_EQ(0.0, intersec[2]);
}

TEST_F(PolygonViewTests, LocalConversionTest) {
    auto query1 = lalib::VecD<3>({ 0.0, 0.5, 0.5 });
    auto local1 = rect->local(query1);
    auto restore1 = rect->point(local1);

    EXPECT_DOUBLE_EQ(restore1[0], query1[0]);
    EXPECT_DOUBLE_EQ(restore1[1], query1[1]);
    EXPECT_DOUBLE_EQ(restore1[2], query1[2]);

    auto query2 = lalib::VecD<3>({ 0.0, 0.2, 0.2 });
    auto local2 = rect2->local(query2);
    auto restore2 = rect2->point(local2);

    EXPECT_DOUBLE_EQ(restore2[0], query2[0]);
    EXPECT_DOUBLE_EQ(restore2[1], query2[1]);
    EXPECT_DOUBLE_EQ(restore2[2], query2[2]);
}

TEST_F(PolygonViewTests, InclusionTest) {
    auto query1 = lalib::VecD<3>({ 0.0, 0.5, 0.5 });
    auto query2 = lalib::VecD<3>({ 0.0, 2.0, 0.5 });

    EXPECT_TRUE(rect->check_inclusion(query1));
    EXPECT_FALSE(rect->check_inclusion(query2));
}
