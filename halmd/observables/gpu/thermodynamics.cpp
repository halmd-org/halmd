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
  // memory allocation in functors
  , sum_velocity_square_()
  , sum_velocity_vector_()
  , sum_scalar_()
  , sum_stress_tensor_diagonal_()
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
double thermodynamics<dimension, float_type>::en_kin()
{
    return .5 * sum_velocity_square_(particle->g_v) / particle->nbox;
}

/**
 * compute mean velocity
 */
template <int dimension, typename float_type>
typename thermodynamics<dimension, float_type>::vector_type
thermodynamics<dimension, float_type>::v_cm()
{
    return sum_velocity_vector_(particle->g_v) / particle->nbox;
}

/**
 * compute potential energy
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::en_pot()
{
    if (!force->aux_flag()) {
        throw std::logic_error("Potential energy not enabled in force module");
    }
    return sum_scalar_(force->potential_energy()) / particle->nbox;
}

/**
 * compute virial sum from potential part of stress tensor
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::virial()
{
    if (!force->aux_flag()) {
        throw std::logic_error("Stress tensor not enabled in force module");
    }
    return sum_stress_tensor_diagonal_(force->stress_tensor_pot()) / particle->nbox;
}

/**
 * compute hypervirial sum
 */
template <int dimension, typename float_type>
double thermodynamics<dimension, float_type>::hypervirial()
{
    if (!force->aux_flag()) {
        throw std::logic_error("Hypervirial not enabled in force module");
    }
    return sum_scalar_(force->hypervirial()) / particle->nbox;
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
