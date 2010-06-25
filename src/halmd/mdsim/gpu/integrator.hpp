/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_INTEGRATOR_HPP
#define HALMD_MDSIM_GPU_INTEGRATOR_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/mdsim/gpu/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <int dimension, typename float_type>
class integrator
  : public mdsim::integrator<dimension>
{
public:
    // module definitions
    typedef integrator _Self;
    typedef mdsim::integrator<dimension> _Base;
    static void options(po::options_description& desc) {}
    static void depends();
    static void select(po::options const& vm) {}

    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef gpu::box<dimension> box_type;
    typedef utility::gpu::device device_type;

    shared_ptr<particle_type> particle;
    shared_ptr<box_type> box;
    shared_ptr<device_type> device;

    integrator(modules::factory& factory, po::options const& vm);
    virtual ~integrator() {}
    virtual void integrate() = 0;
    virtual void finalize() = 0;
};

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATOR_HPP */
