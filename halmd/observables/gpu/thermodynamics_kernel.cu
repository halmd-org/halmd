/*
 * Copyright © 2016      Felix Höfling
 * Copyright © 2013-2015 Nicolas Höft
 * Copyright © 2012      Peter Colberg
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

#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/algorithm/gpu/reduce_kernel.cuh>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/observables/gpu/thermodynamics_kernel.hpp>
#include <halmd/mdsim/force_kernel.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
__device__ void kinetic_energy<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> v;
    float mass;
    tie(v, mass) <<= tex1Dfetch<float4>(texture_, i);
    mv2_ += mass * inner_prod(v, v);
}

template <int dimension, typename float_type>
__device__ void total_force<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> f = tex1Dfetch(force<dimension>::tex_, i);
    force_ += f;
}

template <int dimension, typename float_type>
__device__ void centre_of_mass<dimension, float_type>::operator()(typename iterator::value_type const& value)
{
    size_type i;
    fixed_vector<float, dimension> box_length;
    tie(i, box_length) = value;
    fixed_vector<float, dimension> r, v, img;
    unsigned int species;
    float mass;
    tie(r, species) <<= tex1Dfetch<float4>(position_texture_, i);
    tie(v, mass) <<= tex1Dfetch<float4>(velocity_texture_, i);
    img = tex1Dfetch<coalesced_vector_type>(image_texture_, i);
    mdsim::gpu::box_kernel::extend_periodic(r, img, box_length);
    mr_ += mass * r;
    m_ += mass;
}

template <int dimension, typename float_type>
__device__ void velocity_of_centre_of_mass<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> v;
    float mass;
    tie(v, mass) <<= tex1Dfetch<float4>(texture_, i);
    mv_ += mass * v;
    m_ += mass;
}

template <typename float_type>
__device__ void potential_energy<float_type>::operator()(size_type i)
{
    en_pot_ += tex1Dfetch<float>(texture_, i);
}

template <int dimension, typename float_type>
__device__ void virial<dimension, float_type>::operator()(size_type i)
{
    typedef fixed_vector<float, dimension> stress_pot_diagonal;
    stress_pot_diagonal v;
    v = mdsim::read_stress_tensor_diagonal<stress_pot_diagonal>(texture_, i, stride_);
    // add trace of the potential part of the stress tensor
    for (int j = 0; j < dimension; ++j) {
        virial_ += v[j];
    }
}

template <int dimension, typename float_type>
__device__ void stress_tensor<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> v;
    float mass;

    stress_tensor_ += mdsim::read_stress_tensor<stress_tensor_type>(stress_pot_texture_, i, stride_);
    tie(v, mass) <<= tex1Dfetch<float4>(velocity_texture_, i);
    // compute the kinetic part of the stress tensor
    stress_tensor_ += mass * mdsim::make_stress_tensor(v);
}

template class observables::gpu::kinetic_energy<3, dsfloat>;
template class observables::gpu::kinetic_energy<2, dsfloat>;
template class observables::gpu::total_force<3, dsfloat>;
template class observables::gpu::total_force<2, dsfloat>;
template class observables::gpu::centre_of_mass<3, dsfloat>;
template class observables::gpu::centre_of_mass<2, dsfloat>;
template class observables::gpu::velocity_of_centre_of_mass<3, dsfloat>;
template class observables::gpu::velocity_of_centre_of_mass<2, dsfloat>;
template class observables::gpu::potential_energy<dsfloat>;
template class observables::gpu::virial<3, dsfloat>;
template class observables::gpu::virial<2, dsfloat>;
template class observables::gpu::stress_tensor<3, dsfloat>;
template class observables::gpu::stress_tensor<2, dsfloat>;

} // namespace gpu
} // namespace observables

template class reduction_kernel<observables::gpu::kinetic_energy<3, dsfloat> >;
template class reduction_kernel<observables::gpu::kinetic_energy<2, dsfloat> >;
template class reduction_kernel<observables::gpu::total_force<3, dsfloat> >;
template class reduction_kernel<observables::gpu::total_force<2, dsfloat> >;
template class reduction_kernel<observables::gpu::centre_of_mass<3, dsfloat> >;
template class reduction_kernel<observables::gpu::centre_of_mass<2, dsfloat> >;
template class reduction_kernel<observables::gpu::velocity_of_centre_of_mass<3, dsfloat> >;
template class reduction_kernel<observables::gpu::velocity_of_centre_of_mass<2, dsfloat> >;
template class reduction_kernel<observables::gpu::potential_energy<dsfloat> >;
template class reduction_kernel<observables::gpu::virial<3, dsfloat> >;
template class reduction_kernel<observables::gpu::virial<2, dsfloat> >;
template class reduction_kernel<observables::gpu::stress_tensor<3, dsfloat> >;
template class reduction_kernel<observables::gpu::stress_tensor<2, dsfloat> >;

} // namespace halmd
