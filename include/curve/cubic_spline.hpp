#pragma once 
#include <algorithm>
#include "../../third_party/lalib/include/vec.hpp"
#include "../../third_party/lalib/include/solver/tri_diag.hpp"

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

    auto dp(double s) const noexcept -> VectorType;

    /// @brief  Returns the tangent vector
    /// @param s 
    /// @return 
    auto tangent(double s) const noexcept -> VectorType;

private:
    std::vector<std::array<VectorType, 4>> _c;
    std::vector<double> _cumul_length;

    static auto __calc_coeffs(const std::vector<PointType>& nodes) noexcept 
        -> std::vector<std::array<VectorType, 4>>;

    auto __search_seg_id(double s) const noexcept -> size_t;
    auto __seg_length(size_t seg_id) const noexcept -> double;
    auto __seg_local_pos(double glob_pos) const noexcept -> double;

    auto __calc_seg_p(size_t seg_id, double local_pos) const noexcept -> PointType;
    auto __calc_seg_dp(size_t seg_id, double local_pos) const noexcept -> VectorType;
};

template <size_t N>
inline CubicSpline<N>::CubicSpline(const std::vector<PointType> &nodes):
    _c(__calc_coeffs(nodes)), 
    _cumul_length({0.0})
{
    /// @todo Extract the code section below regarding Simpson's rule 
    // Calculate lengths per segment (Simpson's Rule)
    const auto n_div = 10u;
    auto n_segs = nodes.size() - 1;

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
        this->_cumul_length.emplace_back(total_len);
    }
}

template <size_t N>
inline auto CubicSpline<N>::n_segments() const noexcept -> size_t
{
    auto n = this->_cumul_length.size() - 1;
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
    auto seg_id = this->__search_seg_id(s);
    auto local_s = this->__seg_local_pos(s);

    auto p = this->__calc_seg_p(seg_id, local_s);
    return p;
}

template <size_t N>
inline auto CubicSpline<N>::length() const noexcept -> double
{
    return this->_cumul_length.back();
}

template <size_t N>
inline auto CubicSpline<N>::dp(double s) const noexcept -> VectorType
{
    auto seg_id = this->__search_seg_id(s);
    auto local_s = this->__seg_local_pos(s);

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
inline auto CubicSpline<N>::__search_seg_id(double s) const noexcept -> size_t
{
    auto len = this->length();
    for (auto i = 1u; i < this->_cumul_length.size(); ++i) {
        if (s <= this->_cumul_length[i] / len) {
            return i - 1;
        }
    }
    return this->n_segments() - 1;  // last segment 
}

template <size_t N>
inline auto CubicSpline<N>::__seg_length(size_t seg_id) const noexcept -> double
{
    auto l = this->_cumul_length[seg_id + 1] - this->_cumul_length[seg_id];
    return l;
}

template <size_t N>
inline auto CubicSpline<N>::__seg_local_pos(double glob_pos) const noexcept -> double
{
    auto total_length = this->length();
    auto seg_id = this->__search_seg_id(glob_pos);
    auto s = (glob_pos * total_length - this->_cumul_length[seg_id]) / this->__seg_length(seg_id);
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
}