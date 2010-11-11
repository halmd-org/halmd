/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_GPU_THERMODYNAMICS_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/program_options/program_options.hpp>

namespace halmd
{
namespace observables { namespace gpu
{

template <int dimension, typename float_type>
class thermodynamics
    : public observables::thermodynamics<dimension>
{
public:
    typedef observables::thermodynamics<dimension> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef typename _Base::box_type box_type;
    typedef mdsim::gpu::force<dimension, float_type> force_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<force_type> force;

    static void luaopen(lua_State* L);

    thermodynamics(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<force_type> force
    );

    virtual double en_kin() const;
    virtual vector_type v_cm() const;
    virtual double en_pot() const;
    virtual double virial() const;
};

}} // namespace observables::gpu

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
