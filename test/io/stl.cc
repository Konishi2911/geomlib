#include <gtest/gtest.h>
#include <vector>
#include "geomlib/io/stl.hpp"
#include "geomlib/surface/polygon_mesh.hpp"
#include "lalib/vec.hpp"

using geomlib::Stl;
using geomlib::StlError;
using geomlib::PolygonMesh;

TEST(StlTests, ConstructionTest) {
    auto verts = std::vector {
        lalib::VecD<3>({ 0.0, 0.0, 0.0 }),
        lalib::VecD<3>({ 1.0, 0.0, 0.0 }),
        lalib::VecD<3>({ 1.0, 1.0, 0.0 }),
        lalib::VecD<3>({ 0.0, 1.0, 0.0 }),
        lalib::VecD<3>({ 0.0, 0.0, 1.0 }),
    };
    auto triangular_mesh = PolygonMesh(std::vector(verts), std::vector {
        std::vector<size_t> { 0, 1, 2 },
        std::vector<size_t> { 0, 3, 2 },
        std::vector<size_t> { 0, 4, 2 },
    });
    auto mixed_mesh = PolygonMesh(std::vector(verts), std::vector {
        std::vector<size_t> { 0, 1, 2, 3 },
        std::vector<size_t> { 0, 4, 2 },
    });

    EXPECT_NO_THROW({ 
        auto _ = Stl(std::move(triangular_mesh));
    });
    EXPECT_THROW({ 
        auto _ = Stl(std::move(mixed_mesh));
    }, StlError::IlegalPolygon);
}