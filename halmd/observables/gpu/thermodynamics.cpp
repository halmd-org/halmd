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

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/observables/gpu/thermodynamics.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;
using namespace halmd::algorithm::gpu;

namespace halmd
{
namespace observables { namespace gpu
{

template <int dimension, typename float_type>
thermodynamics<dimension, float_type>::thermodynamics(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<force_type> force
)
  : _Base(box)
  // dependency injection
  , particle(particle)
  , force(force)
{
}

/**
 * preparations before computation of forces
 *
 * set flag of force module to compute auxiliary
 * variables like potential energy, stress tensor,
 * and hypervirial
 */
template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::prepare()
{
    force->aux_enable();
}

/**
 * call sample() from base class and
 * unset flags for auxiliary variables of force module at the end
 */
template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::sample(double time)
{
    _Base::sample(time);
    force->aux_disable();
}

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
    >()(force->potential_energy());

    return en_pot / particle->nbox;
}

/**
 * compute virial sum from potential part of stress tensor
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::virial() const
{
    typedef typename force_type::stress_tensor_type stress_tensor_type;
    typedef typename force_type::gpu_stress_tensor_type gpu_stress_tensor_type;

    // using fixed_vector<dsfloat, N> as output_type results in
    // exit code 255 of 'nvcc -c thermodynamics_kernel.cu'
    double virial = reduce<
        sum_                                      // reduce_transform
      , stress_tensor_type                        // input_type
      , gpu_stress_tensor_type                    // coalesced_input_type
      , dsfloat                                   // output_type
      , dsfloat                                   // coalesced_output_type
      , double                                    // host_output_type
      , at_0                                      // input_transform
    >()(force->stress_tensor_pot());

    return virial / particle->nbox;
}

/**
 * compute hypervirial sum
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::hypervirial() const
{
    double hypervir = reduce<
        sum_                                    // reduce_transform
      , float                                   // input_type
      , float                                   // coalesced_input_type
      , dsfloat                                 // output_type
      , dsfloat                                 // coalesced_output_type
      , double                                  // host_output_type
    >()(force->hypervirial());

    return hypervir / particle->nbox;
}

template <int dimension, typename float_type>
void thermodynamics<dimension, float_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name("thermodynamics_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                namespace_("gpu")
                [
                    class_<thermodynamics, shared_ptr<_Base_Base>, bases<_Base, _Base_Base> >(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                          , shared_ptr<box_type>
                          , shared_ptr<force_type>
                        >())
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
    [
        &thermodynamics<3, float>::luaopen
    ]
    [
        &thermodynamics<2, float>::luaopen
    ];
}

// explicit instantiation
template class thermodynamics<3, float>;
template class thermodynamics<2, float>;

}} // namespace observables::gpu

} // namespace halmd
