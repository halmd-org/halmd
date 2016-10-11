/*
 * Copyright © 2016 Daniel Kirchner
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_SMOOTH_R4_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_SMOOTH_R4_HPP

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/matrix_shape.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * define Lennard-Jones potential and parameters
 */
template <typename potential_type>
class smooth_r4 : public potential_type
{
public:
    typedef typename potential_type::float_type float_type;
    typedef typename potential_type::matrix_type matrix_type;

    template<typename... Args>
    smooth_r4(matrix_type const& cutoff, float_type h, Args&&... args)
            : potential_type (std::forward<Args>(args)...)
            , r_cut_sigma_(check_shape(cutoff, this->sigma()))
            , r_cut_(element_prod(this->sigma(), r_cut_sigma_))
            , rr_cut_(element_prod(r_cut_, r_cut_))
            , en_cut_(this->size1(), this->size2())
            , rri_smooth_(std::pow(h, -2))
    {
        for (size_t i = 0; i < this->size1(); ++i) {
            for (size_t j = 0; j < this->size2(); ++j) {
                std::tie(std::ignore, en_cut_(i,j)) = potential_type::operator()(rr_cut_(i, j), i, j);
            }
        }
        LOG("potential cutoff length: r_c = " << r_cut_sigma_);
        LOG("potential cutoff energy: U = " << en_cut_);
    }

    bool within_range(float_type rr, unsigned a, unsigned b) const
    {
        return rr < rr_cut_(a,b);
    }

    matrix_type const& r_cut() const
    {
        return r_cut_;
    }

    float_type r_cut(unsigned a, unsigned b) const
    {
        return r_cut_(a, b);
    }

    float_type rr_cut(unsigned a, unsigned b) const
    {
        return rr_cut_(a, b);
    }

    matrix_type const& r_cut_sigma() const
    {
        return r_cut_sigma_;
    }

    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        float_type f_abs, pot;
        tie(f_abs, pot) = potential_type::operator()(rr, a, b);
        pot = pot - en_cut_(a,b);
        float_type r = std::sqrt(rr);
        float_type dr = r - this->r_cut(a, b);
        float_type x2 = dr * dr * rri_smooth_;
        float_type x4 = x2 * x2;
        float_type x4i = 1 / (1 + x4);
        // smoothing function
        float_type h0_r = x4 * x4i;
        // first derivative
        float_type h1_r = 4 * dr * rri_smooth_ * x2 * x4i * x4i;
        // apply smoothing function to obtain C¹ force function
        f_abs = h0_r * f_abs - h1_r * (pot / r);
        // apply smoothing function to obtain C² potential function
        pot = h0_r * pot;
        return std::make_tuple(f_abs, pot);
    }
    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L) {
        using namespace luaponte;
        module(L, "libhalmd")
        [
                namespace_("mdsim")
                [
                        namespace_("host")
                        [
                                namespace_("potentials")
                                [
                                        namespace_("pair")
                                        [
                                                class_<smooth_r4, potential_type, std::shared_ptr<smooth_r4> >()
                                                    .property("r_cut", (matrix_type const& (smooth_r4::*)() const) &smooth_r4::r_cut)
                                                    .property("r_cut_sigma", &smooth_r4::r_cut_sigma)
                                              , def("smooth_r4", &std::make_shared<smooth_r4
                                                                , matrix_type const&
                                                                , float_type
                                                                , potential_type const&>)
                                        ]
                                ]
                        ]
                ]
        ];
    }
private:
    /** cutoff length in units of sigma */
    matrix_type r_cut_sigma_;
    /** cutoff length in MD units */
    matrix_type r_cut_;
    /** square of cutoff length */
    matrix_type rr_cut_;
    /** potential energy at cutoff length in MD units */
    matrix_type en_cut_;

    float_type rri_smooth_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_SMOOTH_R4_HPP */
