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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_LOCAL_R4_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_LOCAL_R4_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * define Lennard-Jones potential and parameters
 */
template <typename potential_type>
class local_r4 : public potential_type
{
public:
    typedef typename potential_type::float_type float_type;

    template<typename... Args>
    local_r4(float_type h, Args&&... args)
            : potential_type (std::forward<Args>(args)...), rri_smooth_(std::pow(h, -2)) {
    }

    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        float_type f_abs, pot;
        tie(f_abs, pot) = potential_type::operator()(rr, a, b);
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
                                                class_<local_r4, potential_type, std::shared_ptr<local_r4> >()
                                              , def("local_r4", &std::make_shared<local_r4, float_type, potential_type const&>)
                                        ]
                                ]
                        ]
                ]
        ];
    }
private:
    float_type rri_smooth_;
};

} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_LOCAL_R4_HPP */
