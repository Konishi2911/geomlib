#pragma once 
#ifndef GEOMLIB_CUBIC_SPLINE_HPP
#define GEOMLIB_CUBIC_SPLINE_HPP

#include <algorithm>
#include "../../third_party/lalib/include/vec.hpp"
#include "../../third_party/lalib/include/solver/tri_diag.hpp"
#include "../affine/affine_core.hpp"

namespace geomlib {

template<size_t N>
struct CubicSpline {
public:
    using PointType = lalib::SizedVec<double, N>;
    using VectorType = lalib::SizedVec<double, N>;


    CubicSpline(const std::vector<PointType>& nodes);

    /// @brief  Returns the number of segments in the cubic spline
    /// @return the number of segments
    auto n_segments() const noexcept -> size_t;

    auto operator()(double s) const noexcept -> PointType;

    /// @brief  Returns the point on the segment specified by the given local position `s`
    /// @note   `s` is clamped with 0.0 or 1.0 if it is out of range.
    /// @param s    local position `s` which range is [0.0, 1.0].
    /// @return     a point on the segment
    auto point(double s) const noexcept -> PointType;

    /// @brief  Returns the length of the segment.
    /// @return length of the segment
    auto length() const noexcept -> double;

    /// @brief  Returns the length from starting point to the intermediate position `s`.
    /// @param s 
    auto length(double s) const noexcept -> double;

    auto dp(double s) const noexcept -> VectorType;

    /// @brief  Returns the tangent vector
    /// @param s 
    /// @return 
    auto tangent(double s) const noexcept -> VectorType;

    auto transform(const Affine<N>& affine) noexcept -> CubicSpline<N>&;

    auto transformed(const Affine<N>& affine) const noexcept -> CubicSpline<N>;


private:
    std::vector<std::array<VectorType, 4>> _c;
    std::vector<double> _cumul_length;

    static auto __calc_coeffs(const std::vector<PointType>& nodes) noexcept 
        -> std::vector<std::array<VectorType, 4>>;

    auto __calc_cumul_length(std::vector<double>&) const noexcept -> std::vector<double>&;

    auto __get_seg_id(double s) const noexcept -> size_t;
    auto __seg_length(size_t seg_id) const noexcept -> double;
    auto __seg_local_pos(double glob_pos, size_t segid) const noexcept -> double;

    auto __calc_seg_p(size_t seg_id, double local_pos) const noexcept -> PointType;
    auto __calc_seg_dp(size_t seg_id, double local_pos) const noexcept -> VectorType;

    template<std::invocable<double> F>
    static auto __int(double st, double ed, F&& f, size_t n_div) noexcept -> decltype(f(st));
};

template <size_t N>
inline CubicSpline<N>::CubicSpline(const std::vector<PointType> &nodes):
    _c(__calc_coeffs(nodes)), 
    _cumul_length({0.0})
{
    this->__calc_cumul_length(this->_cumul_length);
}

template <size_t N>
inline auto CubicSpline<N>::n_segments() const noexcept -> size_t
{
    auto n = this->_c.size();
    return n;
}

template <size_t N>
inline auto CubicSpline<N>::operator()(double s) const noexcept -> PointType
{
    return this->point(s);
}

template <size_t N>
inline auto CubicSpline<N>::point(double s) const noexcept -> PointType
{
    auto seg_id = this->__get_seg_id(s);
    auto local_s = this->__seg_local_pos(s, seg_id);

    auto p = this->__calc_seg_p(seg_id, local_s);
    return p;
}

template <size_t N>
inline auto CubicSpline<N>::length() const noexcept -> double
{
    return this->_cumul_length.back();
}

template <size_t N>
inline auto CubicSpline<N>::length(double s) const noexcept -> double
{
    auto sid = this->__get_seg_id(s);
    auto seg_s = this->__seg_local_pos(s, sid);
    auto len = this->_cumul_length[sid];
    len += __int(0.0, seg_s, [&](double s){ return this->__calc_seg_dp(sid, s).norm2(); }, 10);
    return len;
}

template <size_t N>
inline auto CubicSpline<N>::dp(double s) const noexcept -> VectorType
{
    auto seg_id = this->__get_seg_id(s);
    auto local_s = this->__seg_local_pos(s, seg_id);

    auto p = this->__calc_seg_dp(seg_id, local_s);
    return p;
}

template <size_t N>
inline auto CubicSpline<N>::tangent(double s) const noexcept -> VectorType
{
    auto p = this->dp(s);
    p = p / p.norm2();
    return p;
}

template <size_t N>
inline auto CubicSpline<N>::transform(const Affine<N> &affine) noexcept -> CubicSpline<N> &
{
    auto n = this->_c.size();
    for (auto i = 0u; i < n; ++i) {
        affine.transform(this->_c[i][0]);
        affine.transform(this->_c[i][1]);
        affine.transform(this->_c[i][2]);
        affine.transform(this->_c[i][3]);
    }
    this->__calc_cumul_length(this->_cumul_length);
    return *this;
}

template <size_t N>
inline auto CubicSpline<N>::transformed(const Affine<N> &affine) const noexcept -> CubicSpline<N>
{
    auto curve = CubicSpline(*this).transform(affine);
    return curve;
}

// ==== Private Functions ==== //

template <size_t N>
inline auto CubicSpline<N>::__calc_coeffs(const std::vector<PointType> &nodes) noexcept 
    -> std::vector<std::array<VectorType, 4>>
{
    // Create coefficient matrix and RHS matrix
    auto n = nodes.size() - 1;
    auto mat = lalib::DynTriDiagMat<double>(
        std::vector<double>(n - 2, 1.0),
        std::vector<double>(n - 1, 4.0),
        std::vector<double>(n - 2, 1.0)
    );
    auto tmp_rhs = std::vector<double>();
    tmp_rhs.reserve((n - 1) * N);
    for (auto i = 0u; i < n - 1; ++i) {
        auto p = 3.0 * (nodes[i + 2] - 2.0 * nodes[i + 1] + nodes[i]);
        tmp_rhs.insert(tmp_rhs.end(), p.begin(), p.end());
    }
    auto rhs = lalib::DynMat<double>(std::move(tmp_rhs), n - 1, N);

    // Solve linear equation of Tridiagonal Matrix 
    auto solver = lalib::solver::TriDiag(std::move(mat));
    solver.solve_linear(rhs, rhs);

    /// @todo   Consider changing to `vector-view` to avoid unnecessary copy.
    auto c = std::vector<lalib::VecD<N>>();
    c.reserve(n - 1);
    for (auto i = 0u; i < n - 1; ++i) {
        auto tmp = std::array<double, N>();
        for (auto j = 0u; j < N; ++j) {
            tmp[j] = rhs(i, j);
        }
        c.emplace_back(tmp);
    }

    // Create segment-by-segment coefficient.
    auto coeffs = std::vector<std::array<VectorType, 4>>(); 
    coeffs.emplace_back(
        std::array{
            nodes[0],
            nodes[1] - nodes[0] - c[0] / 3.0,
            lalib::SizedVec<double, N>::filled(0.0),
            c[0] / 3.0
        }
    );
    for (auto i = 1u; i < n - 1; ++i) {
        coeffs.emplace_back(
            std::array{
                nodes[i],
                nodes[i + 1] - nodes[i] - (c[i] + 2.0 * c[i - 1]) / 3.0,
                c[i - 1],
                (c[i] - c[i - 1]) / 3.0
            }
        );
    }
    coeffs.emplace_back(
        std::array{
            nodes[n - 1],
            nodes[n] - nodes[n - 1] - 2.0 * c[n - 2] / 3.0,
            c[n - 2],
            - c[n - 2] / 3.0
        }
    );

    return coeffs;
}

template <size_t N>
inline auto CubicSpline<N>::__calc_cumul_length(std::vector<double>& cumul_length) const noexcept -> std::vector<double>&
{
    /// @todo Extract the code section below regarding Simpson's rule 
    // Calculate lengths per segment (Simpson's Rule)
    const auto n_div = 10u;
    auto n_segs = this->_c.size();

    double total_len = 0.0;
    for (auto i = 0u; i < n_segs; ++i) {
        double len = 0.0;
        for (auto k = 0u; k <= n_div - 2; k += 2) {
            const auto dp0 = this->__calc_seg_dp(i, static_cast<double>(k) / n_div).norm2();
            const auto dp1 = this->__calc_seg_dp(i, static_cast<double>(k + 1) / n_div).norm2();
            const auto dp2 = this->__calc_seg_dp(i, static_cast<double>(k + 2) / n_div).norm2();

            len += (dp0 + 4.0 * dp1 + dp2) / (3.0 * n_div);
        }
        total_len += len;
        cumul_length.emplace_back(total_len);
    }
    return cumul_length;
}

template <size_t N>
inline auto CubicSpline<N>::__get_seg_id(double s) const noexcept -> size_t
{
    s = std::clamp(s, 0.0, 1.0 - std::numeric_limits<double>::epsilon());
    auto n_segs = this->n_segments();
    auto sid = static_cast<size_t>(std::floor(s * n_segs));
    return sid;
}

template <size_t N>
inline auto CubicSpline<N>::__seg_length(size_t seg_id) const noexcept -> double
{
    auto l = this->_cumul_length[seg_id + 1] - this->_cumul_length[seg_id];
    return l;
}

template <size_t N>
inline auto CubicSpline<N>::__seg_local_pos(double glob_pos, size_t segid) const noexcept -> double
{
    auto n_segs = this->n_segments();
    auto s = glob_pos * n_segs - segid;
    return s;
}

template <size_t N>
inline auto CubicSpline<N>::__calc_seg_p(size_t seg_id, double local_s) const noexcept -> PointType
{
    auto p = this->_c[seg_id][0] + this->_c[seg_id][1] * local_s + this->_c[seg_id][2] * std::pow(local_s, 2) + this->_c[seg_id][3] * std::pow(local_s, 3);
    return p;
}

template <size_t N>
inline auto CubicSpline<N>::__calc_seg_dp(size_t seg_id, double s) const noexcept -> VectorType
{
    auto p = this->_c[seg_id][1] + 2.0 * this->_c[seg_id][2] * s + 3.0 * this->_c[seg_id][3] * std::pow(s, 2);
    return p;
}

template <size_t N>
template <std::invocable<double> F>
inline auto CubicSpline<N>::__int(double st, double ed, F &&f, size_t n_div) noexcept -> decltype(f(st))
{
    auto ds = ed - st;
    auto len = decltype(f(std::declval<double>()))();
    for (auto k = 0u; k <= n_div - 2; k += 2) {
        const auto s0 = ds * static_cast<double>(k) / n_div + st;
        const auto s1 = ds * static_cast<double>(k + 1) / n_div + st;
        const auto s2 = ds * static_cast<double>(k + 2) / n_div + st;
        const auto f0 = f(s0);
        const auto f1 = f(s1);
        const auto f2 = f(s2);

        len += (f0 + 4.0 * f1 + f2) * ds / (3.0 * n_div);
    }
    return len;
}

}

#endif