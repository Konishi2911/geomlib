#pragma once

#include "geomlib/surface/polygon.hpp"
#include "geomlib/surface/surface_concepts.hpp"

#include <numeric>
#include <span>

namespace geomlib {

struct PolygonMesh {
    using PointType = lalib::VecD<3>;
    using VectorType = lalib::VecD<3>;

    PolygonMesh(std::vector<lalib::VecD<3>>&& verts, const std::vector<std::vector<size_t>>& vlinks) noexcept;

    auto area() const noexcept -> double;
    auto n_faces() const noexcept -> size_t;
    auto verts() const noexcept -> std::span<const PointType>;
    auto subfaces() const noexcept -> std::span<const PolygonView>;

private:
    std::vector<lalib::VecD<3>> _verts;
    std::vector<PolygonView> _elems;

    double _area;

    static auto __create_elements(const std::vector<lalib::VecD<3>>& verts, const std::vector<std::vector<size_t>>& vlinks) -> std::vector<PolygonView>;
    static auto __calc_area(const std::vector<PolygonView>&) noexcept -> double;
};


inline PolygonMesh::PolygonMesh(std::vector<lalib::VecD<3>>&& verts, const std::vector<std::vector<size_t>>& vlinks) noexcept:
    _verts(std::move(verts)),
    _elems(__create_elements(_verts, vlinks)),
    _area(__calc_area(_elems))
{}


inline auto PolygonMesh::area() const noexcept -> double {
    return this->_area;
}

inline auto PolygonMesh::n_faces() const noexcept -> size_t {
    return this->_elems.size();
}

inline auto PolygonMesh::verts() const noexcept -> std::span<const PointType> {
    return std::span(this->_verts);
}

inline auto PolygonMesh::subfaces() const noexcept -> std::span<const PolygonView> {
    return std::span(this->_elems);
}


inline auto PolygonMesh::__create_elements(const std::vector<lalib::VecD<3>>& verts, const std::vector<std::vector<size_t>>& vlinks) -> std::vector<PolygonView> {
    auto elems = std::vector<PolygonView>();
    elems.reserve(vlinks.size());
    for (auto i = 0u; i < vlinks.size(); ++i) {
        auto tmp = std::vector<std::reference_wrapper<const lalib::VecD<3>>>();
        for (auto&& vid: vlinks[i]) {
            tmp.emplace_back(std::cref(verts[vid]));
        }
        elems.emplace_back(std::move(tmp));
    }
    return elems; 
}

inline auto PolygonMesh::__calc_area(const std::vector<PolygonView>& faces) noexcept -> double {
    auto area = std::accumulate(faces.begin(), faces.end(), 0.0, [](auto&& acc, auto&& face){ return acc + face.area(); });
    return area;
}

}
