#include <gtest/gtest.h>
#include <numbers>
#include "geomlib/affine.hpp"
#include "geomlib/curve/segment.hpp"

TEST(AffineTests, Rotation2DTest) {
    const auto affine = geomlib::rotate2d(std::numbers::pi / 6.0, lalib::VecD<2>::filled(0.0));
    auto ex = lalib::VecD<2>({ 1.0, 0.0 });
    auto ey = lalib::VecD<2>({ 0.0, 1.0 });
    
    affine.transform(ex, ex);
    affine.transform(ey);

    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ex[0]);
    ASSERT_DOUBLE_EQ(0.5 , ex[1]);

    ASSERT_DOUBLE_EQ(-0.5 , ey[0]);
    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ey[1]);
}

TEST(AffineTests, Rotation2DWithPivotTest) {
    auto pivot = lalib::VecD<2>({ 1.0, 0.0 });
    const auto affine = geomlib::rotate2d(std::numbers::pi / 6.0, pivot);
    auto zero = lalib::VecD<2>({ 0.0, 0.0 });
    auto ex = lalib::VecD<2>({ 1.0, 0.0 });

    affine.transform(zero, zero);
    ex = affine.transformed(ex);

    ASSERT_DOUBLE_EQ(1.0 - std::sqrt(3) / 2.0, zero[0]);
    ASSERT_DOUBLE_EQ(-0.5, zero[1]);

    ASSERT_DOUBLE_EQ(1.0, ex[0]);
    ASSERT_DOUBLE_EQ(0.0, ex[1]);
}

TEST(AffineTests, Rotation3DXTest) {
    const auto affine = geomlib::rotate<3, 0>(std::numbers::pi / 6.0);
    auto ex = lalib::VecD<3>({ 1.0, 0.0, 0.0 });
    auto ey = lalib::VecD<3>({ 0.0, 1.0, 0.0 });
    auto ez = lalib::VecD<3>({ 0.0, 0.0, 1.0 });
    
    affine.transform(ex, ex);
    affine.transform(ey, ey);
    affine.transform(ez, ez);

    ASSERT_DOUBLE_EQ(1.0, ex[0]);
    ASSERT_DOUBLE_EQ(0.0, ex[1]);
    ASSERT_DOUBLE_EQ(0.0, ex[2]);

    ASSERT_DOUBLE_EQ(0.0 , ey[0]);
    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ey[1]);
    ASSERT_DOUBLE_EQ(0.5 , ey[2]);

    ASSERT_DOUBLE_EQ(0.0 , ez[0]);
    ASSERT_DOUBLE_EQ(-0.5 , ez[1]);
    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ez[2]);
}

TEST(AffineTests, Rotation3DYTest) {
    const auto affine = geomlib::rotate<3, 1>(std::numbers::pi / 6.0);
    auto ex = lalib::VecD<3>({ 1.0, 0.0, 0.0 });
    auto ey = lalib::VecD<3>({ 0.0, 1.0, 0.0 });
    auto ez = lalib::VecD<3>({ 0.0, 0.0, 1.0 });
    
    affine.transform(ex, ex);
    affine.transform(ey, ey);
    affine.transform(ez, ez);

    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ex[0]);
    ASSERT_DOUBLE_EQ(0.0 , ex[1]);
    ASSERT_DOUBLE_EQ(-0.5 , ex[2]);

    ASSERT_DOUBLE_EQ(0.0, ey[0]);
    ASSERT_DOUBLE_EQ(1.0, ey[1]);
    ASSERT_DOUBLE_EQ(0.0, ey[2]);

    ASSERT_DOUBLE_EQ(0.5 , ez[0]);
    ASSERT_DOUBLE_EQ(0.0 , ez[1]);
    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ez[2]);
}

TEST(AffineTests, Rotation3DZTest) {
    const auto affine = geomlib::rotate<3, 2>(std::numbers::pi / 6.0);
    auto ex = lalib::VecD<3>({ 1.0, 0.0, 0.0 });
    auto ey = lalib::VecD<3>({ 0.0, 1.0, 0.0 });
    auto ez = lalib::VecD<3>({ 0.0, 0.0, 1.0 });
    
    affine.transform(ex, ex);
    affine.transform(ey, ey);
    affine.transform(ez, ez);

    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ex[0]);
    ASSERT_DOUBLE_EQ(0.5 , ex[1]);
    ASSERT_DOUBLE_EQ(0.0 , ex[2]);

    ASSERT_DOUBLE_EQ(-0.5 , ey[0]);
    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ey[1]);
    ASSERT_DOUBLE_EQ(0.0 , ey[2]);

    ASSERT_DOUBLE_EQ(0.0, ez[0]);
    ASSERT_DOUBLE_EQ(0.0, ez[1]);
    ASSERT_DOUBLE_EQ(1.0, ez[2]);
}

TEST(AffineTEsts, Rotation3DTest) {
    auto ex = lalib::VecD<3>({ 1.0, 0.0, 0.0 });
    auto ey = lalib::VecD<3>({ 0.0, 1.0, 0.0 });
    auto ez = lalib::VecD<3>({ 0.0, 0.0, 1.0 });
    const auto affine = geomlib::rotate3d(std::numbers::pi / 6.0, ez, lalib::VecD<3>::filled(0.0));
    
    affine.transform(ex, ex);
    affine.transform(ey, ey);
    affine.transform(ez, ez);

    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ex[0]);
    ASSERT_DOUBLE_EQ(0.5 , ex[1]);
    ASSERT_DOUBLE_EQ(0.0 , ex[2]);

    ASSERT_DOUBLE_EQ(-0.5 , ey[0]);
    ASSERT_DOUBLE_EQ(std::sqrt(3) / 2.0 , ey[1]);
    ASSERT_DOUBLE_EQ(0.0 , ey[2]);

    ASSERT_DOUBLE_EQ(0.0, ez[0]);
    ASSERT_DOUBLE_EQ(0.0, ez[1]);
    ASSERT_DOUBLE_EQ(1.0, ez[2]);
}

TEST(AffineTests, SegmentRotationTest) {
    auto seg = geomlib::Segment(
        lalib::VecD<2>({ 0.0, 0.0 }),
        lalib::VecD<2>({ 1.0, 1.0 })
    );
    auto affine = geomlib::rotate2d(std::numbers::pi / 4.0, lalib::VecD<2>::filled(0.0));

    auto trans_seg = geomlib::transform_lazy(std::move(seg), std::move(affine));
    
    ASSERT_DOUBLE_EQ(0.0, trans_seg.point(0.0)[0]);
    ASSERT_DOUBLE_EQ(0.0, trans_seg.point(0.0)[1]);

    ASSERT_NEAR(0.0, trans_seg.point(1.0)[0], 1.0e-14);
    ASSERT_DOUBLE_EQ(std::sqrt(2.0), trans_seg.point(1.0)[1]);
}

TEST(AffineTests, CompositeTest) {
    auto p = lalib::VecD<2>({2.0, 0.0});
    auto a1 = geomlib::rotate2d(std::numbers::pi / 4.0, lalib::VecD<2>({ 1.0, 0.0 }));
    auto a2 = geomlib::rotate2d(std::numbers::pi / 4.0, lalib::VecD<2>({ 1.0, 0.0 }));

    a1.composite(a2);
    a1.transform(p);

    ASSERT_DOUBLE_EQ(1.0, p[0]);
    ASSERT_DOUBLE_EQ(1.0, p[1]);
}


TEST(AffineTransformedCurveTests, MultipleLazyTransformationTest) {
    auto curve = geomlib::Segment<2>(
        lalib::VecD<2>({0.0, 0.0}),
        lalib::VecD<2>({2.0, 0.0})
    );
    auto affine = geomlib::rotate2d(std::numbers::pi / 4.0, lalib::VecD<2>::filled(0.0));

	static_assert(geomlib::AffineTransformableCurve<decltype(curve), 2>);

    auto transformed_curve = geomlib::transform_lazy(std::move(curve), std::move(affine));
    transformed_curve = geomlib::transform_lazy(std::move(transformed_curve), std::move(affine));

    EXPECT_DOUBLE_EQ(0.0, transformed_curve.point(0.0)[0]);
    EXPECT_DOUBLE_EQ(0.0, transformed_curve.point(0.0)[1]);

    EXPECT_NEAR(0.0, transformed_curve.point(1.0)[0], 1e-10);
    EXPECT_NEAR(2.0, transformed_curve.point(1.0)[1], 1e-10); 
}