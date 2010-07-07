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

#ifndef HALMD_MDSIM_GPU_INTEGRATORS_VERLET_HPP
#define HALMD_MDSIM_GPU_INTEGRATORS_VERLET_HPP

#include <cuda_wrapper.hpp>

#include <halmd/mdsim/gpu/integrator.hpp>
#include <halmd/mdsim/gpu/integrators/verlet_kernel.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/util/exception.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

template <int dimension, typename float_type>
class verlet
  : public mdsim::gpu::integrator<dimension, float_type>
{
public:
    // module definitions
    typedef verlet _Self;
    typedef mdsim::gpu::integrator<dimension, float_type> _Base;
    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::options const& vm);

    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::particle_type particle_type;
    typedef typename _Base::box_type box_type;

    using _Base::particle;
    using _Base::box;
    using _Base::device;

    /** CUDA C++ wrapper */
    verlet_wrapper<dimension> const* wrapper;

    verlet(modules::factory& factory, po::options const& vm);
    virtual ~verlet() {}
    void integrate();
    void finalize();

protected:
    /** integration time-step */
    using _Base::timestep_;
    /** half time-step */
    double timestep_half_;
};

}}} // namespace mdsim::gpu::integrators

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATORS_VERLET_HPP */
