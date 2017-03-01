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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_SHARP_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_SHARP_HPP

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
namespace truncations {

/**
 * define sharp potential truncation
 */
template <typename potential_type>
class sharp
  : public potential_type
{
public:
    typedef typename potential_type::float_type float_type;
    typedef typename potential_type::matrix_type matrix_type;

    template<typename... Args>
    sharp(matrix_type const& cutoff, Args&&... args)
            : potential_type (std::forward<Args>(args)...)
            , r_cut_sigma_(check_shape(cutoff, this->sigma()))
            , r_cut_(element_prod(this->sigma(), r_cut_sigma_))
            , rr_cut_(element_prod(r_cut_, r_cut_))
    {
        LOG("potential cutoff length: r_c = " << r_cut_sigma_);
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
        return potential_type::operator()(rr, a, b);
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L)
    {
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
                            class_<sharp, potential_type, std::shared_ptr<sharp> >()
                                .property("r_cut", (matrix_type const& (sharp::*)() const) &sharp::r_cut)
                                .property("r_cut_sigma", &sharp::r_cut_sigma)

                          , def("sharp", &std::make_shared<sharp
                              , matrix_type const&
                              , potential_type const&>
                            )
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
};

} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_SHARP_HPP */
