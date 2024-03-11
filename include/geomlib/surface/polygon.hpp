#pragma once
#include "geomlib/surface/surface_concepts.hpp"
#ifndef GEOMLIB_SURFACE_POLYGON_HPP
#define GEOMLIB_SURFACE_POLYGON_HPP

#include "lalib/vec.hpp"
#include "lalib/mat.hpp"

namespace geomlib {

struct Polygon {
    using PointType = lalib::VecD<3>;
    using VectorType = lalib::VecD<3>;

    Polygon(std::vector<PointType>&& verts) noexcept;

    auto point(const lalib::VecD<2>& s) const -> PointType;
    auto normal(const lalib::VecD<2>& s = lalib::VecD<2>::filled(0.0)) const -> VectorType;

    auto verts() const noexcept -> const std::vector<PointType>&;
    auto origin() const noexcept -> PointType;
    auto area() const noexcept -> double;
    auto base() const noexcept -> lalib::MatD<3, 2>;

private:
    std::vector<PointType> _verts;
    std::array<lalib::VecD<3>, 2> _base;

    static auto __calc_base(const std::vector<PointType>& verts) noexcept -> std::array<lalib::VecD<3>, 2>;
};


inline Polygon::Polygon(std::vector<PointType>&& verts) noexcept: 
    _verts(verts), _base(__calc_base(verts))
{}

inline auto Polygon::verts() const noexcept -> const std::vector<PointType>& {
    return this->_verts;
}

inline auto Polygon::point(const lalib::VecD<2>& s) const -> PointType {
    auto p = this->origin();
    lalib::axpy(s[0], this->_base[0], p);
    lalib::axpy(s[1], this->_base[1], p);
    return p;
}

inline auto Polygon::normal(const lalib::VecD<2>&) const -> PointType {
    auto n = lalib::cross(this->_base[0], this->_base[1]);
    return n;
}

inline auto Polygon::origin() const noexcept -> PointType {
    auto o = this->_verts[0];
    return o;
}

inline auto Polygon::area() const noexcept -> double {
    auto darea = lalib::VecD<3>::filled(0.0);
    for (auto i = 1u; i < this->_verts.size() - 1; ++i) {
        auto p1 = this->_verts[i] - this->_verts[0];
        auto p2= this->_verts[i + 1] - this->_verts[0];
        lalib::cross_acc(p1, p2, darea);
    }
    auto area = 0.5 * darea.norm2();
    return area;
}

inline auto Polygon::base() const noexcept -> lalib::MatD<3, 2> {
    auto base = lalib::MatD<3, 2>({
        this->_base[0][0], this->_base[1][0],
        this->_base[0][1], this->_base[1][1],
        this->_base[0][2], this->_base[1][2]
    });
    return base;
}


inline auto Polygon::__calc_base(const std::vector<PointType>& verts) noexcept -> std::array<lalib::VecD<3>, 2> {
    auto n = verts.size();
    auto u1 = verts[1] - verts[0]; 
    auto u2 = verts[n - 1] - verts[0]; 
    u1 = u1 / u1.norm2();
    u2 = u2 / u2.norm2();

    return std::array { u1, u2 };
}

static_assert(Surface<Polygon>);

}

#endif