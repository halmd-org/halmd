/*
 * Copyright © 2014 Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_EXTERNAL_HARMONIC_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_EXTERNAL_HARMONIC_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <lua.hpp>
#include <memory>
#include <tuple>

#include <halmd/io/logger.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace external {

/**
 * define harmonic potential and parameters
 */
template <int dimension, typename float_type>
class harmonic
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef boost::numeric::ublas::vector<float_type> scalar_container_type;
    typedef boost::numeric::ublas::vector<vector_type> vector_container_type;

    harmonic(
        scalar_container_type const& stiffness
      , vector_container_type const& offset
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Compute force and potential for interaction.
     *
     * @param r   particle position reduced to periodic box
     * @param species   particle species
     * @returns   tuple of force vector @f$ -\nabla U(\vec r) @f$ and potential @f$ U(\vec r) @f$
     */
    std::tuple<vector_type, float_type> operator()(vector_type const& r, unsigned int species) const
    {
        vector_type dr = r - offset_[species];
        vector_type force = - stiffness_[species] * dr;
        float_type en_pot = - inner_prod(force, dr) / 2;

        return std::make_tuple(force, en_pot);
    }

    scalar_container_type const& stiffness() const
    {
        return stiffness_;
    }

    vector_container_type const& offset() const
    {
        return offset_;
    }

    unsigned int size() const
    {
        return stiffness_.size();
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** potential stiffness in MD units */
    scalar_container_type stiffness_;
    /** position of potential minimum in MD units */
    vector_container_type offset_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace external
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_EXTERNAL_HARMONIC_HPP */
