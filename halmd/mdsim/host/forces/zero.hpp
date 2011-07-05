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

#ifndef HALMD_MDSIM_HOST_FORCES_ZERO_HPP
#define HALMD_MDSIM_HOST_FORCES_ZERO_HPP

#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

/**
 * zero force (noninteracting particles)
 */
template <int dimension, typename float_type>
class zero
  : public mdsim::host::force<dimension, float_type>
{
public:
    typedef mdsim::host::force<dimension, float_type> _Base;
    typedef typename _Base::stress_tensor_type stress_tensor_type;

    typedef host::particle<dimension, float_type> particle_type;

    static char const* module_name() { return "zero"; }

    boost::shared_ptr<particle_type> particle;

    static void luaopen(lua_State* L);

    zero(boost::shared_ptr<particle_type> particle);

    // there's nothing to compute
    virtual void compute() {}

    // nothing to enable or disable
    virtual void aux_enable() {}
    virtual void aux_disable() {}
    virtual bool aux_flag() const
    {
        return true;
    }

    //! return average potential energy per particle
    virtual double potential_energy()
    {
        return 0;
    }

    //! potential part of stress tensor
    virtual stress_tensor_type stress_tensor_pot()
    {
        return stress_tensor_type(0);
    }

    //! return hypervirial per particle
    virtual double hypervirial()
    {
        return 0;
    }
};

} // namespace mdsim
} // namespace host
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_ZERO_HPP */
