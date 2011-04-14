/*
 * Copyright © 2010-2011  Felix Höfling
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

#ifndef HALMD_OBSERVABLES_HOST_THERMODYNAMICS_HPP
#define HALMD_OBSERVABLES_HOST_THERMODYNAMICS_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/observables/thermodynamics.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>

namespace halmd
{
namespace observables { namespace host
{

template <int dimension, typename float_type>
class thermodynamics
    : public observables::thermodynamics<dimension>
{
public:
    typedef observables::thermodynamics<dimension> _Base;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef typename _Base::box_type box_type;
    typedef typename _Base::core_type core_type;
    typedef mdsim::host::force<dimension, float_type> force_type;
    typedef typename particle_type::vector_type vector_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<force_type> force;

    static void luaopen(lua_State* L);

    thermodynamics(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<core_type> core
      , boost::shared_ptr<force_type> force
    );

    virtual void prepare();
    virtual void sample(uint64_t step);

    virtual double en_kin();
    virtual vector_type v_cm();

    virtual double en_pot()
    {
        if (!force->aux_flag()) {
            throw std::logic_error("Potential energy not enabled in force module");
        }
        return force->potential_energy();
    }

    virtual double virial()
    {
        if (!force->aux_flag()) {
            throw std::logic_error("Stress tensor not enabled in force module");
        }
        return force->stress_tensor_pot()[0];
    }

    virtual double hypervirial()
    {
        if (!force->aux_flag()) {
            throw std::logic_error("Hypervirial not enabled in force module");
        }
        return force->hypervirial();
    }
};

}} // namespace observables::host

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_THERMODYNAMICS_HPP */
