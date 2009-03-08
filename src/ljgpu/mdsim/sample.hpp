/* Phase space sample
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef LJGPU_MDSIM_SAMPLE_HPP
#define LJGPU_MDSIM_SAMPLE_HPP

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#ifdef WITH_CUDA
# include <cuda_wrapper.hpp>
#endif
#include <ljgpu/mdsim/impl.hpp>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/math/vector4d.hpp>
#include <vector>

namespace ljgpu {

/**
 * phase space sample for binary mixture
 */
template <typename float_type, int dimension>
struct trajectory_host_sample
{
    typedef vector<double, dimension> position_vector;
    typedef vector<float_type, dimension> velocity_vector;
    typedef std::vector<position_vector> position_sample_vector;
    typedef std::vector<velocity_vector> velocity_sample_vector;
    typedef boost::shared_ptr<position_sample_vector> position_sample_ptr;
    typedef boost::shared_ptr<velocity_sample_vector> velocity_sample_ptr;

    trajectory_host_sample() {}
    trajectory_host_sample(position_sample_ptr r, velocity_sample_ptr v) : r(r), v(v) {}

    /** periodically extended particle positions */
    position_sample_ptr r;
    /** particle velocities */
    velocity_sample_ptr v;
};

#if WITH_CUDA

template <int dimension>
struct trajectory_gpu_sample;

template <>
struct trajectory_gpu_sample<3>
{
    typedef float4 position_vector;
    typedef float4 velocity_vector;
    typedef cuda::vector<position_vector> position_sample_vector;
    typedef cuda::vector<velocity_vector> velocity_sample_vector;
    typedef boost::shared_ptr<position_sample_vector> position_sample_ptr;
    typedef boost::shared_ptr<velocity_sample_vector> velocity_sample_ptr;

    trajectory_gpu_sample() {}
    trajectory_gpu_sample(position_sample_ptr r, velocity_sample_ptr v) : r(r), v(v) {}

    /** periodically extended particle positions */
    position_sample_ptr r;
    /** particle velocities */
    velocity_sample_ptr v;
};

template <>
struct trajectory_gpu_sample<2>
{
    typedef float2 position_vector;
    typedef float2 velocity_vector;
    typedef cuda::vector<position_vector> position_sample_vector;
    typedef cuda::vector<velocity_vector> velocity_sample_vector;
    typedef boost::shared_ptr<position_sample_vector> position_sample_ptr;
    typedef boost::shared_ptr<velocity_sample_vector> velocity_sample_ptr;

    trajectory_gpu_sample() {}
    trajectory_gpu_sample(position_sample_ptr r, velocity_sample_ptr v) : r(r), v(v) {}

    /** periodically extended particle positions */
    position_sample_ptr r;
    /** particle velocities */
    velocity_sample_ptr v;
};

#endif /* WITH_CUDA */

/**
 * MD simulation sample for A and B particles
 */
template <int dimension>
struct energy_sample
{
    typedef vector<double, 1 + (dimension - 1) * dimension / 2> virial_tensor;

    /** mean potential energy per particle */
    double en_pot;
    /** virial tensor trace and off-diagonal elements for particle species */
    std::vector<virial_tensor> virial;
    /** mean squared velocity per particle */
    double vv;
    /** mean velocity per particle */
    vector<double, dimension> v_cm;
};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_SAMPLE_HPP */
