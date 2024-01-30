#include "../../include/curve/segment.hpp"
#include <gtest/gtest.h>

TEST(SegmentTests, LengthTest) {
    auto sp = lalib::VecD<3>({0.0, 1.0, 2.0});
    auto ep = lalib::VecD<3>({2.0, 4.0, 1.0});
    auto seg = geomlib::Segment(sp, ep);

    ASSERT_DOUBLE_EQ(std::sqrt(14.0), seg.length());
}

TEST(SegmentTests, InterPointTest) {
    auto sp = lalib::VecD<3>({0.0, 1.0, 2.0});
    auto ep = lalib::VecD<3>({2.0, 4.0, 1.0});
    auto seg = geomlib::Segment(sp, ep);

    // Edge Cases
    EXPECT_DOUBLE_EQ(sp[0], seg.point(0.0)[0]);
    EXPECT_DOUBLE_EQ(sp[1], seg.point(0.0)[1]);
    EXPECT_DOUBLE_EQ(sp[2], seg.point(0.0)[2]);

    EXPECT_DOUBLE_EQ(ep[0], seg.point(1.0)[0]);
    EXPECT_DOUBLE_EQ(ep[1], seg.point(1.0)[1]);
    EXPECT_DOUBLE_EQ(ep[2], seg.point(1.0)[2]);

    // Internal points
    EXPECT_DOUBLE_EQ(1.0, seg.point(0.5)[0]);
    EXPECT_DOUBLE_EQ(2.5, seg.point(0.5)[1]);
    EXPECT_DOUBLE_EQ(1.5, seg.point(0.5)[2]);

    // Out-of-range points
    EXPECT_DOUBLE_EQ(ep[0], seg.point(2.0)[0]);
    EXPECT_DOUBLE_EQ(ep[1], seg.point(2.0)[1]);
    EXPECT_DOUBLE_EQ(ep[2], seg.point(2.0)[2]);
}

TEST(SegmentTests, TangentTest) {
    auto sp = lalib::VecD<3>({0.0, 1.0, 2.0});
    auto ep = lalib::VecD<3>({2.0, 4.0, 1.0});
    auto seg = geomlib::Segment(sp, ep);

    EXPECT_DOUBLE_EQ(2.0/std::sqrt(14.0), seg.tangent(0.0)[0]);
    EXPECT_DOUBLE_EQ(3.0/std::sqrt(14.0), seg.tangent(0.0)[1]);
    EXPECT_DOUBLE_EQ(-1.0/std::sqrt(14.0), seg.tangent(0.0)[2]);
}

TEST(SegmentTests, DistanceTest) {
    auto sp = lalib::VecD<2>({0.0, 0.0});
    auto ep = lalib::VecD<2>({1.0, 0.0});
    auto seg = geomlib::Segment(sp, ep);

    auto q1 = lalib::VecD<2>({0.5, 1.0});
    auto q2 = lalib::VecD<2>({-1.0, 0.0});
    auto q3 = lalib::VecD<2>({2.0, 1.0});

    EXPECT_DOUBLE_EQ(1.0, seg.distance(q1));
    EXPECT_DOUBLE_EQ(1.0, seg.distance(q2));
    EXPECT_DOUBLE_EQ(std::sqrt(2.0), seg.distance(q3));
}

TEST(SegmentTests, AffineTransformationTest) {
    auto sp = lalib::VecD<3>({0.0, 1.0, 2.0});
    auto ep = lalib::VecD<3>({2.0, 4.0, 1.0});
    auto seg = geomlib::Segment(sp, ep);

    auto affine = geomlib::Affine<3>(
        lalib::MatD<3, 3>({
            1.0, 0.0, 0.0,
            0.0, 0.5, 0.0,
            0.0, 0.0, 1.0
        }),
        lalib::VecD<3>::filled(0.0)
    );

    seg.transform(affine);

    EXPECT_DOUBLE_EQ(std::sqrt(7.25), seg.length());
}