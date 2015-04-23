/*
 * Copyright © 2013-2015 Nicolas Höft
 * Copyright © 2012      Peter Colberg
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

#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/algorithm/gpu/reduce_kernel.cuh>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/observables/gpu/thermodynamics_kernel.hpp>
#include <halmd/mdsim/force_kernel.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/** particle positions and species */
static texture<float4> position_;
/** particle velocities and masses */
static texture<float4> velocity_;
/** potential energies */
static texture<float> en_pot_;
/** potential parts of stress tensor */
static texture<float> stress_pot_;

/** minimum image vectors */
template<int dimension>
struct image
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> image<dimension>::type image<dimension>::tex_;

template <int dimension, typename float_type>
__device__ void kinetic_energy<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> v;
    float mass;
    tie(v, mass) <<= tex1Dfetch(velocity_, i);
    mv2_ += mass * inner_prod(v, v);
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
    tie(r, species) <<= tex1Dfetch(position_, i);
    tie(v, mass) <<= tex1Dfetch(velocity_, i);
    img = tex1Dfetch(image<dimension>::tex_, i);
    mdsim::gpu::box_kernel::extend_periodic(r, img, box_length);
    mr_ += mass * r;
    m_ += mass;
}

template <int dimension, typename float_type>
__device__ void velocity_of_centre_of_mass<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> v;
    float mass;
    tie(v, mass) <<= tex1Dfetch(velocity_, i);
    mv_ += mass * v;
    m_ += mass;
}

template <typename float_type>
__device__ void potential_energy<float_type>::operator()(size_type i)
{
    en_pot_ += tex1Dfetch(gpu::en_pot_, i);
}

template <int dimension, typename float_type>
__device__ void virial<dimension, float_type>::operator()(size_type i)
{
    typedef fixed_vector<float, dimension> stress_pot_diagonal;
    stress_pot_diagonal v;
    v = mdsim::read_stress_tensor_diagonal<stress_pot_diagonal>(stress_pot_, i, stride_);
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

    stress_tensor_ += mdsim::read_stress_tensor<stress_tensor_type>(stress_pot_, i, stride_);
    tie(v, mass) <<= tex1Dfetch(velocity_, i);
    // compute the kinetic part of the stress tensor
    stress_tensor_ += mass * mdsim::make_stress_tensor(v);
}

template <int dimension, typename float_type>
cuda::texture<float4> const
kinetic_energy<dimension, float_type>::texture_ = velocity_;

template <int dimension, typename float_type>
cuda::texture<float4> const
centre_of_mass<dimension, float_type>::position_texture_ = position_;

template <int dimension, typename float_type>
cuda::texture<typename centre_of_mass<dimension, float_type>::coalesced_vector_type> const
centre_of_mass<dimension, float_type>::image_texture_ = image<dimension>::tex_;

template <int dimension, typename float_type>
cuda::texture<float4> const
centre_of_mass<dimension, float_type>::velocity_texture_ = velocity_;

template <int dimension, typename float_type>
cuda::texture<float4> const
velocity_of_centre_of_mass<dimension, float_type>::texture_ = velocity_;

template <typename float_type>
cuda::texture<float> const
potential_energy<float_type>::texture_ = gpu::en_pot_;

template <int dimension, typename float_type>
cuda::texture<float> const
virial<dimension, float_type>::texture_ = stress_pot_;

template <int dimension, typename float_type>
cuda::texture<float> const
stress_tensor<dimension, float_type>::stress_pot_texture_ = stress_pot_;

template <int dimension, typename float_type>
cuda::texture<float4> const
stress_tensor<dimension, float_type>::velocity_texture_ = velocity_;

template class observables::gpu::kinetic_energy<3, dsfloat>;
template class observables::gpu::kinetic_energy<2, dsfloat>;
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
