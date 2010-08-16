/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/mdsim/gpu/thermodynamics.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;
using namespace std;
using namespace halmd::algorithm::gpu;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::depends()
{
    modules::depends<_Self, particle_type>::required();
}

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(factory, vm))
  // allocate result variables
  , virial_(particle->ntype)
  , g_en_pot_(particle->dim.threads())
  , g_virial_(particle->dim.threads())
{}

/**
 * compute mean-square velocity
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_kin() const
{
    double vv = reduce<
        sum_                                    // reduce_transform
      , fixed_vector<float, dimension>          // input_type
      , float4                                  // coalesced_input_type
      , dsfloat                                 // output_type
      , dsfloat                                 // coalesced_output_type
      , double                                  // host_output_type
      , square_                                 // input_transform
    >()(particle->g_v);

    return .5 * vv / particle->nbox;
}

/**
 * compute mean velocity
 */
template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type
thermodynamics<dimension, float_type>::v_cm() const
{
    vector_type v = reduce<
        sum_                                    // reduce_transform
      , fixed_vector<float, dimension>          // input_type
      , float4                                  // coalesced_input_type
      , fixed_vector<dsfloat, dimension>        // output_type
      , fixed_vector<dsfloat, dimension>        // coalesced_output_type
      , vector_type                             // host_output_type
    >()(particle->g_v);

    return v / particle->nbox;
}

/**
 * compute potential energy
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_pot() const
{
    double en_pot = reduce<
        sum_                                    // reduce_transform
      , float                                   // input_type
      , float                                   // coalesced_input_type
      , dsfloat                                 // output_type
      , dsfloat                                 // coalesced_output_type
      , double                                  // host_output_type
    >()(g_en_pot_);

    return en_pot / particle->nbox;
}

/**
 * compute virial sum regardless of particle type
 */
template <int dimension, typename float_type>
std::vector<typename thermodynamics<dimension, float_type>::virial_type> const&
thermodynamics<dimension, float_type>::virial()
{
    typedef fixed_vector<float, virial_type::static_size> float_virial_type;

    // using fixed_vector<dsfloat, N> as output_type results in
    // exit code 255 of 'nvcc -c thermodynamics_kernel.cu'
    virial_[0] = reduce<
        sum_                               // reduce_transform
      , float_virial_type                  // input_type
      , gpu_virial_type                    // coalesced_input_type
      , float_virial_type                  // output_type
      , float_virial_type                  // coalesced_output_type
      , virial_type                        // host_output_type
    >()(g_virial_);

    virial_[0] /= particle->nbox;

    return virial_;
}

// explicit instantiation
template class thermodynamics<3, float>;
template class thermodynamics<2, float>;

}} // namespace mdsim::gpu

template class module<mdsim::gpu::thermodynamics<3, float> >;
template class module<mdsim::gpu::thermodynamics<2, float> >;

} // namespace halmd
