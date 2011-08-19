/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_FORCES_ZERO_HPP
#define HALMD_MDSIM_GPU_FORCES_ZERO_HPP

#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/particle.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * zero force (noninteracting particles)
 */
template <int dimension, typename float_type>
class zero
  : public mdsim::gpu::force<dimension, float_type>
{
public:
    typedef mdsim::gpu::force<dimension, float_type> _Base;
    typedef typename _Base::gpu_stress_tensor_type gpu_stress_tensor_type;

    typedef gpu::particle<dimension, float_type> particle_type;

    static char const* module_name() { return "zero"; }

    boost::shared_ptr<particle_type> particle;

    static void luaopen(lua_State* L);
    /** zero potential energy, stress tensor and hypervirial */
    zero(boost::shared_ptr<particle_type> particle);
    /** zero particle forces */
    virtual void compute();

    // nothing to enable or disable
    virtual void aux_enable() {}
    virtual void aux_disable() {}
    virtual bool aux_flag() const
    {
        return true;
    }

    //! returns potential energies of particles
    virtual cuda::vector<float> const& potential_energy() const
    {
        return g_en_pot_;
    }

    /** potential part of stress tensors of particles */
    virtual cuda::vector<gpu_stress_tensor_type> const& stress_tensor_pot() const
    {
        return g_stress_pot_;
    }

    //! returns hyper virial of particles
    virtual cuda::vector<float> const& hypervirial() const
    {
        return g_hypervirial_;
    }

protected:
    /** potential energy for each particle */
    cuda::vector<float> g_en_pot_;
    /** potential part of stress tensor for each particle */
    cuda::vector<gpu_stress_tensor_type> g_stress_pot_;
    /** hyper virial for each particle */
    cuda::vector<float> g_hypervirial_;
};

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_ZERO_HPP */
