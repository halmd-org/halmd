/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_hoover.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

template <int dimension, typename float_type>
verlet_nvt_hoover<dimension, float_type>::
verlet_nvt_hoover(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , float_type timestep
  , float_type temperature
  , float_type resonance_frequency
)
  // dependency injection
  : particle(particle)
  , box(box)
  // member initialisation
  , xi(0)
  , v_xi(0)
  , en_nhc_(0)
  , resonance_frequency_(resonance_frequency)
  // set up functor and allocate internal memory,
  // this is done here and only once rather than repeatedly during the integration
  , compute_en_kin_2_(30)  // FIXME reduce kernel should use #multiprocessors = #blocks as default
{
    this->timestep(timestep);

    LOG("resonance frequency of heat bath: " << resonance_frequency_);
    this->temperature(temperature);

#ifdef USE_VERLET_DSFUN
    //
    // Double-single precision requires two single precision
    // "words" per coordinate. We use the first part of a GPU
    // vector for the higher (most significant) words of all
    // particle positions or velocities, and the second part for
    // the lower (least significant) words.
    //
    // The additional memory is allocated using reserve(), which
    // increases the capacity() without changing the size().
    //
    // Take care to pass capacity() as an argument to cuda::copy
    // or cuda::memset calls if needed, as the lower words will
    // be ignored in the operation.
    //
    LOG("using velocity-Verlet integration in double-single precision");
    particle->g_r.reserve(2 * particle->dim.threads());
    // particle images remain in single precision as they
    // contain integer values (and otherwise would not matter
    // for the long-time stability of the Verlet integrator)
    particle->g_v.reserve(2 * particle->dim.threads());
#else
    LOG_WARNING("using velocity-Verlet integration in single precision");
#endif

    // copy parameters to CUDA device
    try {
        cuda::copy(
            static_cast<fixed_vector<float, dimension> >(box->length())
          , wrapper_type::kernel.box_length
        );
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to initialize Verlet integrator symbols");
        throw;
    }
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.integrate, "integrate", "first half-step of velocity-Verlet (+ Nosé-Hoover chain)");
    profiler.register_runtime(runtime_.finalize, "finalize", "second half-step of velocity-Verlet (+ Nosé-Hoover chain)");
    profiler.register_runtime(runtime_.propagate, "propagate", "propagate Nosé-Hoover chain");
    profiler.register_runtime(runtime_.rescale, "rescale", "rescale velocities in Nosé-Hoover thermostat");
}

/**
 * register observables
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::register_observables(writer_type& writer)
{
    writer.register_observable("ENHC", &en_nhc_, "energy of Nosé-Hoover chain variables per particle");
}

/**
 * set integration time-step
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
timestep(double timestep)
{
    timestep_ = static_cast<float_type>(timestep);
    timestep_half_ = timestep_ / 2;
    timestep_4_ = timestep_ / 4;
    timestep_8_ = timestep_ / 8;

    try {
        cuda::copy(timestep_, wrapper_type::kernel.timestep);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to initialize Verlet integrator symbols");
        throw;
    }

    LOG("integration timestep: " << timestep_);
}

/*
 * set temperature and adjust masses of heat bath variables
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
temperature(double temperature)
{
    temperature_ = static_cast<float_type>(temperature);
    en_kin_target_2_ = dimension * particle->nbox * temperature_;

    // follow Martyna et al. [J. Chem. Phys. 97, 2635 (1992)]
    // for the masses of the heat bath variables
    float_type omega_sq = pow(2 * M_PI * resonance_frequency_, 2);
    unsigned int dof = dimension * particle->nbox;
    fixed_vector<float_type, 2> mass;
    mass[0] = dof * temperature_ / omega_sq;
    mass[1] = temperature_ / omega_sq;
    this->mass(mass);

    LOG("temperature of heat bath: " << temperature_);
    LOG_DEBUG("target kinetic energy: " << en_kin_target_2_ / particle->nbox);
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
mass(fixed_vector<double, 2> const& mass)
{
    mass_xi_ = static_cast<fixed_vector<float_type, 2> >(mass);
    LOG("`mass' of heat bath variables: " << mass_xi_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
integrate()
{
    scoped_timer<timer> timer_(runtime_.integrate);
    float_type scale = propagate_chain();

    try {
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper_type::kernel.integrate(
            particle->g_r, particle->g_image, particle->g_v, particle->g_f, scale
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream first leapfrog step on GPU");
        throw;
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::
finalize()
{
    scoped_timer<timer> timer_(runtime_.finalize);

    // TODO: possibly a performance critical issue:
    // the old implementation had this loop included in update_forces(),
    // which saves one additional read of the forces plus the additional kernel execution
    // and scheduling
    try {
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper_type::kernel.finalize(particle->g_v, particle->g_f);
        cuda::thread::synchronize();

        float_type scale = propagate_chain();

        // rescale velocities
        scoped_timer<timer> timer2_(runtime_.rescale);
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper_type::kernel.rescale(particle->g_v, scale);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream second leapfrog step on GPU");
        throw;
    }

    // compute energy contribution of chain variables
    en_nhc_ = temperature_ * (dimension * particle->nbox * xi[0] + xi[1]);
    for (unsigned int i = 0; i < 2; ++i ) {
        en_nhc_ += mass_xi_[i] * v_xi[i] * v_xi[i] / 2;
    }
    en_nhc_ /= particle->nbox;
}

/**
 * propagate Nosé-Hoover chain
 */
template <int dimension, typename float_type>
float_type verlet_nvt_hoover<dimension, float_type>::propagate_chain()
{
    scoped_timer<timer> timer_(runtime_.propagate);

    // compute total kinetic energy (multiplied by 2),
    float_type en_kin_2 = compute_en_kin_2_(particle->g_v);

    // head of the chain
    v_xi[1] += (mass_xi_[0] * v_xi[0] * v_xi[0] - temperature_) / mass_xi_[1] * timestep_4_;
    float_type t = exp(-v_xi[1] * timestep_8_);
    v_xi[0] *= t;
    v_xi[0] += (en_kin_2 - en_kin_target_2_) / mass_xi_[0] * timestep_4_;
    v_xi[0] *= t;

    // propagate heat bath variables
    for (unsigned int i = 0; i < 2; ++i ) {
        xi[i] += v_xi[i] * timestep_half_;
    }

    // rescale velocities and kinetic energy
    // we only compute the factor, the rescaling is done elsewhere
    float_type s = exp(-v_xi[0] * timestep_half_);
    en_kin_2 *= s * s;

    // tail of the chain, mirrors the head
    v_xi[0] *= t;
    v_xi[0] += (en_kin_2 - en_kin_target_2_) / mass_xi_[0] * timestep_4_;
    v_xi[0] *= t;
    v_xi[1] += (mass_xi_[0] * v_xi[0] * v_xi[0] - temperature_) / mass_xi_[1] * timestep_4_;

    // return scaling factor for CUDA kernels
    return s;
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(verlet_nvt_hoover<dimension, float_type> const&)
{
    return verlet_nvt_hoover<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void verlet_nvt_hoover<dimension, float_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("integrators")
                [
                    class_<
                        verlet_nvt_hoover
                      , shared_ptr<_Base_Base>
                      , bases<_Base_Base, _Base>
                    >(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type>
                          , shared_ptr<box_type>
                          , float_type, float_type, float_type
                        >())
                        .def("register_runtimes", &verlet_nvt_hoover::register_runtimes)
                        .def("register_observables", &verlet_nvt_hoover::register_observables)
                        .property("mass", (fixed_vector<double, 2> const& (verlet_nvt_hoover::*)() const)&verlet_nvt_hoover::mass) // FIXME make read/write
                        .property("resonance_frequency", &verlet_nvt_hoover::resonance_frequency)
                        .property("en_nhc", &verlet_nvt_hoover::en_nhc)
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_verlet_nvt_hoover(lua_State* L)
{
#ifdef USE_VERLET_DSFUN
    verlet_nvt_hoover<3, double>::luaopen(L);
    verlet_nvt_hoover<2, double>::luaopen(L);
#else
    verlet_nvt_hoover<3, float>::luaopen(L);
    verlet_nvt_hoover<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_VERLET_DSFUN
template class verlet_nvt_hoover<3, double>;
template class verlet_nvt_hoover<2, double>;
#else
template class verlet_nvt_hoover<3, float>;
template class verlet_nvt_hoover<2, float>;
#endif

}}} // namespace mdsim::gpu::integrators

} // namespace halmd
