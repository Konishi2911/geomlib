#include <gtest/gtest.h>
#include "geomlib/surface/polygon_mesh.hpp"

using geomlib::PolygonMesh;

class PolygonMeshTests: public ::testing::Test {
protected:
    std::vector<lalib::VecD<3>> points = std::vector({
        lalib::VecD<3>({ 0.0, 0.0, 0.0 }),
        lalib::VecD<3>({ 1.0, 0.0, 0.0 }),
        lalib::VecD<3>({ 0.0, 1.0, 0.0 }),
        lalib::VecD<3>({ 0.0, 0.0, 1.0 }),
        lalib::VecD<3>({ 1.0, 1.0, 0.0 }),
        lalib::VecD<3>({ 0.0, 1.0, 1.0 }),
        lalib::VecD<3>({ 1.0, 1.0, 1.0 }),
        lalib::VecD<3>({ 1.0, 0.0, 1.0 }),
    });

    virtual void SetUp() override {

    }
};

TEST_F(PolygonMeshTests, NElementsTest) {
    auto vlinks = std::vector {
        std::vector<size_t> { 0, 1, 4, 2 },
        std::vector<size_t> { 0, 1, 4, 2 },
        std::vector<size_t> { 0, 1, 4, 2 }
    };
    auto mesh = PolygonMesh(std::vector(points), vlinks);
    
    ASSERT_EQ(vlinks.size(), mesh.n_faces());
}

TEST_F(PolygonMeshTests, ElementTest) {
    auto vlinks = std::vector {
        std::vector<size_t> { 0, 1, 4, 2 },
        std::vector<size_t> { 0, 1, 3 },
        std::vector<size_t> { 3, 7, 1 }
    };
    auto mesh = PolygonMesh(std::vector(points), vlinks);
    
    ASSERT_EQ(vlinks.size(), mesh.n_faces());

    auto e0 = mesh.subfaces()[0];
    ASSERT_DOUBLE_EQ(1.0, e0.area());

    auto e1 = mesh.subfaces()[1];
    ASSERT_DOUBLE_EQ(0.5, e1.area());

    auto e2 = mesh.subfaces()[2];
    ASSERT_DOUBLE_EQ(0.5, e2.area());

    ASSERT_DOUBLE_EQ(e0.area() + e1.area() + e2.area(), mesh.area());
}