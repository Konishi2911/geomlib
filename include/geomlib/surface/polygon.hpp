#pragma once
#ifndef GEOMLIB_SURFACE_POLYGON_HPP
#define GEOMLIB_SURFACE_POLYGON_HPP

#include "geomlib/surface/surface_concepts.hpp"
#include "lalib/ops/vec_ops.hpp"
#include "lalib/vec.hpp"
#include "lalib/mat.hpp"
#include <functional>

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

    /// @brief Calculates a foot of a perpendicular
    auto foot(const lalib::VecD<3>& query) const noexcept -> lalib::VecD<3>;

    /// @brief calculates the position in the local coordinate from the given point
    auto local(const VectorType& query) const noexcept -> lalib::VecD<2>;

    /// @brief Checks inclusion of the point
    auto check_inclusion(const lalib::VecD<3>& query) const noexcept -> bool;

private:
    std::vector<PointType> _verts;
    std::array<lalib::VecD<3>, 2> _base;

    static auto __calc_base(const std::vector<PointType>& verts) noexcept -> std::array<lalib::VecD<3>, 2>;
};


struct PolygonView {
    using PointType = lalib::VecD<3>;
    using VectorType = lalib::VecD<3>;

    PolygonView(std::vector<std::reference_wrapper<const PointType>>&& verts) noexcept;

    auto point(const lalib::VecD<2>& s) const -> PointType;
    auto normal(const lalib::VecD<2>& s = lalib::VecD<2>::filled(0.0)) const -> VectorType;

    auto verts() const noexcept -> const std::vector<std::reference_wrapper<const PointType>>&;
    auto origin() const noexcept -> PointType;
    auto area() const noexcept -> double;
    auto base() const noexcept -> lalib::MatD<3, 2>;

    /// @brief Calculates a foot of a perpendicular
    auto foot(const lalib::VecD<3>& query) const noexcept -> lalib::VecD<3>;

    /// @brief calculates the position in the local coordinate from the given point
    auto local(const VectorType& query) const noexcept -> lalib::VecD<2>;

    /// @brief  Checks if the query point is inside of the polygon
    /// @note   This function DO NOT check whether the query point is on the plane belonging to this polygon.
    auto check_inclusion(const lalib::VecD<3>& query) const noexcept -> bool;

private:
    std::vector<std::reference_wrapper<const PointType>> _verts;
    std::array<lalib::VecD<3>, 2> _base;

    static auto __calc_base(const std::vector<std::reference_wrapper<const PointType>>& verts) noexcept -> std::array<lalib::VecD<3>, 2>;
};


// =========================== //
// Polygon                     //
// =========================== //

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
    n = n / n.norm2();
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

inline auto Polygon::foot(const lalib::VecD<3>& query) const noexcept -> lalib::VecD<3> {
    auto pf = query - (query - this->_verts[0]).dot(this->normal()) * this->normal();
    return pf;
}

inline auto Polygon::local(const VectorType& query) const noexcept -> lalib::VecD<2> {
    auto t_arr = std::array<double, 6>();
    for (auto i = 0u; i < 3; ++i) {
        t_arr[i] = this->_base[0][i];
        t_arr[i + 3] = this->_base[1][i];
    }
    auto t = lalib::MatD<2, 3>(std::move(t_arr));
    auto uv = this->_base[0].dot(this->_base[1]);
    auto tt_inv = lalib::MatD<2, 2>({
        1.0 , -uv, -uv, 1.0
    });
    auto tt = lalib::MatD<2, 3>::uninit();
    lalib::mul(1.0 / (1.0 - std::pow(uv, 2)), tt_inv, t, 1.0, tt);
    
    auto q = tt * query;
    return q;
}

inline auto Polygon::check_inclusion(const lalib::VecD<3>& query) const noexcept -> bool {
    auto l0 = (this->_verts[1] - this->_verts[0]);
    auto cross0 = lalib::cross(query - this->_verts[0], l0);
    auto n = this->_verts.size();

    bool check = true;
    for (auto i = 1u; i < n; ++i) {
        auto l = this->_verts[(i + 1) % n] - this->_verts[i];
        auto cross = lalib::cross(query - this->_verts[i], l);
        check &= cross0.dot(cross) > 0;
    }
    return check;
}


inline auto Polygon::__calc_base(const std::vector<PointType>& verts) noexcept -> std::array<lalib::VecD<3>, 2> {
    auto n = verts.size();
    auto u1 = verts[1] - verts[0]; 
    auto u2 = verts[n - 1] - verts[0]; 
    u1 = u1 / u1.norm2();
    u2 = u2 / u2.norm2();

    return std::array { u1, u2 };
}

// =========================== //
// PolygonView                 //
// =========================== //

inline PolygonView::PolygonView(std::vector<std::reference_wrapper<const PointType>>&& verts) noexcept: 
    _verts(verts), _base(__calc_base(verts))
{}

inline auto PolygonView::verts() const noexcept -> const std::vector<std::reference_wrapper<const PointType>>& {
    return this->_verts;
}

inline auto PolygonView::point(const lalib::VecD<2>& s) const -> PointType {
    auto p = this->origin();
    lalib::axpy(s[0], this->_base[0], p);
    lalib::axpy(s[1], this->_base[1], p);
    return p;
}

inline auto PolygonView::normal(const lalib::VecD<2>&) const -> PointType {
    auto n = lalib::cross(this->_base[0], this->_base[1]);
    n = n / n.norm2();
    return n;
}

inline auto PolygonView::origin() const noexcept -> PointType {
    auto o = this->_verts[0];
    return o;
}

inline auto PolygonView::area() const noexcept -> double {
    auto darea = lalib::VecD<3>::filled(0.0);
    for (auto i = 1u; i < this->_verts.size() - 1; ++i) {
        auto p1 = this->_verts[i].get() - this->_verts[0].get();
        auto p2= this->_verts[i + 1].get() - this->_verts[0].get();
        lalib::cross_acc(p1, p2, darea);
    }
    auto area = 0.5 * darea.norm2();
    return area;
}

inline auto PolygonView::base() const noexcept -> lalib::MatD<3, 2> {
    auto base = lalib::MatD<3, 2>({
        this->_base[0][0], this->_base[1][0],
        this->_base[0][1], this->_base[1][1],
        this->_base[0][2], this->_base[1][2]
    });
    return base;
}

inline auto PolygonView::foot(const lalib::VecD<3>& query) const noexcept -> lalib::VecD<3> {
    auto pf = query - (query - this->_verts[0].get()).dot(this->normal()) * this->normal();
    return pf;
}

inline auto PolygonView::local(const VectorType& query) const noexcept -> lalib::VecD<2> {
    auto t_arr = std::array<double, 6>();
    for (auto i = 0u; i < 3; ++i) {
        t_arr[i] = this->_base[0][i];
        t_arr[i + 3] = this->_base[1][i];
    }
    auto t = lalib::MatD<2, 3>(std::move(t_arr));
    auto uv = this->_base[0].dot(this->_base[1]);
    auto tt_inv = lalib::MatD<2, 2>({
        1.0 , -uv, -uv, 1.0
    });
    auto tt = lalib::MatD<2, 3>::uninit();
    lalib::mul(1.0 / (1.0 - std::pow(uv, 2)), tt_inv, t, 1.0, tt);
    
    auto q = tt * query;
    return q;
}

inline auto PolygonView::check_inclusion(const lalib::VecD<3>& query) const noexcept -> bool {
    auto l0 = (this->_verts[1].get() - this->_verts[0].get());
    auto cross0 = lalib::cross(query - this->_verts[0].get(), l0);
    auto n = this->_verts.size();

    bool check = true;
    for (auto i = 1u; i < n; ++i) {
        auto l = this->_verts[(i + 1) % n].get() - this->_verts[i].get();
        auto cross = lalib::cross(query - this->_verts[i].get(), l);
        check &= cross0.dot(cross) > 0;
    }
    return check;
}


inline auto PolygonView::__calc_base(const std::vector<std::reference_wrapper<const PointType>>& verts) noexcept -> std::array<lalib::VecD<3>, 2> {
    auto n = verts.size();
    auto u1 = verts[1].get() - verts[0].get();
    auto u2 = verts[n - 1].get() - verts[0].get();
    u1 = u1 / u1.norm2();
    u2 = u2 / u2.norm2();

    return std::array { u1, u2 };
}

static_assert(Surface<Polygon>);
static_assert(Surface<PolygonView>);

}

#endif