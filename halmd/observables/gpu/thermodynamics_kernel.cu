/*
 * Copyright © 2012  Peter Colberg
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

#include <halmd/algorithm/gpu/reduce_kernel.cuh>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/observables/gpu/thermodynamics_kernel.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/** particle velocities and masses */
static texture<float4> velocity_;
/** potential energies */
static texture<float> en_pot_;
/** potential parts of stress tensor */
static texture<void> stress_pot_;

template <int dimension, typename float_type>
void kinetic_energy<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> v;
    float mass;
    tie(v, mass) <<= tex1Dfetch(velocity_, i);
    mv2_ += mass * inner_prod(v, v);
}

template <int dimension, typename float_type>
void velocity_of_centre_of_mass<dimension, float_type>::operator()(size_type i)
{
    fixed_vector<float, dimension> v;
    float mass;
    tie(v, mass) <<= tex1Dfetch(velocity_, i);
    mv_ += mass * v;
    m_ += mass;
}

template <typename float_type>
void potential_energy<float_type>::operator()(size_type i)
{
    en_pot_ += tex1Dfetch(gpu::en_pot_, i);
}

template <int dimension, typename float_type>
void virial<dimension, float_type>::operator()(size_type i)
{
    stress_pot_type s = tex1Dfetch(reinterpret_cast<texture<coalesced_stress_pot_type>&>(stress_pot_), i);
    virial_ += s[0];
}

template <int dimension, typename float_type>
cuda::texture<float4> const
kinetic_energy<dimension, float_type>::texture_ = velocity_;

template <int dimension, typename float_type>
cuda::texture<float4> const
velocity_of_centre_of_mass<dimension, float_type>::texture_ = velocity_;

template <typename float_type>
cuda::texture<float> const
potential_energy<float_type>::texture_ = gpu::en_pot_;

template <int dimension, typename float_type>
cuda::texture<typename virial<dimension, float_type>::coalesced_stress_pot_type> const
virial<dimension, float_type>::texture_ = stress_pot_;

template class observables::gpu::kinetic_energy<3, dsfloat>;
template class observables::gpu::kinetic_energy<2, dsfloat>;
template class observables::gpu::velocity_of_centre_of_mass<3, dsfloat>;
template class observables::gpu::velocity_of_centre_of_mass<2, dsfloat>;
template class observables::gpu::potential_energy<dsfloat>;
template class observables::gpu::virial<3, dsfloat>;
template class observables::gpu::virial<2, dsfloat>;

} // namespace gpu
} // namespace observables

template class reduction_kernel<observables::gpu::kinetic_energy<3, dsfloat> >;
template class reduction_kernel<observables::gpu::kinetic_energy<2, dsfloat> >;
template class reduction_kernel<observables::gpu::velocity_of_centre_of_mass<3, dsfloat> >;
template class reduction_kernel<observables::gpu::velocity_of_centre_of_mass<2, dsfloat> >;
template class reduction_kernel<observables::gpu::potential_energy<dsfloat> >;
template class reduction_kernel<observables::gpu::virial<3, dsfloat> >;
template class reduction_kernel<observables::gpu::virial<2, dsfloat> >;

} // namespace halmd
