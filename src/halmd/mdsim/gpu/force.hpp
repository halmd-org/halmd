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

#ifndef HALMD_MDSIM_GPU_FORCE_HPP
#define HALMD_MDSIM_GPU_FORCE_HPP

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/gpu/box.hpp>
// #include <halmd/mdsim/gpu/forces/smooth.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace observables { namespace gpu
{

// forward declaration
template <int dimension, typename float_type>
class thermodynamics;

}} // namespace observables::gpu

namespace mdsim { namespace gpu
{

template <int dimension, typename float_type>
class force
  : public mdsim::force<dimension>
{
public:
    // module definitions
    typedef force _Self;
    typedef mdsim::force<dimension> _Base;
    static void options(po::options_description& desc) {}
    static void depends();
    static void select(po::options const& vm) {}

    typedef fixed_vector<float_type, dimension> vector_type;
    typedef boost::numeric::ublas::symmetric_matrix<float_type, boost::numeric::ublas::lower> matrix_type;
    typedef fixed_vector<float_type, 1 + (dimension - 1) * dimension / 2> stress_tensor_type;
    typedef typename boost::mpl::if_c<
        stress_tensor_type::static_size >= 3
      , float4
      , float2
    >::type gpu_stress_tensor_type;

    typedef gpu::particle<dimension, float> particle_type;
    typedef gpu::box<dimension> box_type;
//     typedef host::forces::smooth<dimension, float_type> smooth_type;

    shared_ptr<particle_type> particle;
    shared_ptr<box_type> box;
//     shared_ptr<smooth_type> smooth;

    force(modules::factory& factory, po::options const& vm);
    virtual ~force() {};
    virtual void compute() = 0;
    virtual matrix_type const& cutoff() = 0;

protected:
    /** potential energy for each particle */
    cuda::vector<float> g_en_pot_;
    /** potential part of stress tensor for each particle */
    cuda::vector<gpu_stress_tensor_type> g_stress_pot_;

    friend class observables::gpu::thermodynamics<dimension, float_type>;
};

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCE_HPP */
