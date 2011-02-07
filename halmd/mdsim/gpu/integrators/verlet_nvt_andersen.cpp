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

#include <algorithm>
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using boost::fusion::at_key;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

template <int dimension, typename float_type, typename RandomNumberGenerator>
verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
verlet_nvt_andersen(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<random_type> random
  , float_type timestep, float_type temperature, float_type coll_rate
)
  // dependency injection
  : particle(particle)
  , box(box)
  , random(random)
  , coll_rate_(coll_rate)
{
    this->timestep(timestep);
    this->temperature(temperature);
    LOG("collision rate with heat bath: " << coll_rate_);

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
        cuda::copy(static_cast<vector_type>(box->length()), wrapper_type::kernel.box_length);
        cuda::copy(random->rng.rng(), wrapper_type::kernel.rng);
    }
    catch (cuda::error const& e) {
        LOG_ERROR(e.what());
        throw runtime_error("failed to initialize Verlet integrator symbols");
    }
}

/**
 * set integration time-step
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
timestep(double timestep)
{
    timestep_ = timestep;
    timestep_half_ = 0.5 * timestep_;
    coll_prob_ = coll_rate_ * timestep;

    try {
        cuda::copy(timestep_, wrapper_type::kernel.timestep);
        cuda::copy(coll_prob_, wrapper_type::kernel.coll_prob);
    }
    catch (cuda::error const& e) {
        LOG_ERROR(e.what());
        throw runtime_error("failed to initialize Verlet integrator symbols");
    }

    LOG("integration timestep: " << timestep_);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
temperature(double temperature)
{
  temperature_ = static_cast<float_type>(temperature);
  sqrt_temperature_ = sqrt(temperature_);

  try {
      cuda::copy(sqrt_temperature_, wrapper_type::kernel.sqrt_temperature);
  }
  catch (cuda::error const& e) {
      LOG_ERROR(e.what());
      throw runtime_error("failed to initialize Verlet integrator symbols");
  }

  LOG("temperature of heat bath: " << temperature_);
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
integrate()
{
    try {
        scoped_timer<timer> timer_(at_key<integrate_>(runtime_));
        cuda::configure(particle->dim.grid, particle->dim.block);
        wrapper_type::kernel.integrate(
            particle->g_r, particle->g_image, particle->g_v, particle->g_f);
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::runtime_error("failed to stream first leapfrog step on GPU");
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
finalize()
{
    // TODO: possibly a performance critical issue:
    // the old implementation had this loop included in update_forces(),
    // which saves one additional read of the forces plus the additional kernel execution
    // and scheduling
    try {
        scoped_timer<timer> timer_(at_key<finalize_>(runtime_));
        // use CUDA execution dimensions of 'random' since
        // the kernel makes use of the random number generator
        cuda::configure(random->rng.dim.grid, random->rng.dim.block);
        wrapper_type::kernel.finalize(
            particle->g_v, particle->g_f
          , particle->nbox, particle->dim.threads()
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw std::runtime_error("failed to stream second leapfrog step on GPU");
    }
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
static char const* module_name_wrapper(verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator> const&)
{
    return verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::module_name();
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void verlet_nvt_andersen<dimension, float_type, RandomNumberGenerator>::
luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    namespace_("integrators")
                    [
                        class_<
                            verlet_nvt_andersen
                          , shared_ptr<_Base_Base>
                          , bases<_Base_Base, _Base>
                        >(class_name.c_str())
                            .def(constructor<
                                shared_ptr<particle_type>
                              , shared_ptr<box_type>
                              , shared_ptr<random_type>
                              , float_type, float_type, float_type>()
                            )
                            .def("register_runtimes", &verlet_nvt_andersen::register_runtimes)
                            .property("collision_rate", &verlet_nvt_andersen::collision_rate)
                            .property("module_name", &module_name_wrapper<dimension, float_type, RandomNumberGenerator>)
                    ]
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
    [
        &verlet_nvt_andersen<3, float, random::gpu::rand48>::luaopen
    ]
    [
        &verlet_nvt_andersen<2, float, random::gpu::rand48>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class verlet_nvt_andersen<3, float, random::gpu::rand48>;
template class verlet_nvt_andersen<2, float, random::gpu::rand48>;

}}} // namespace mdsim::gpu::integrators

} // namespace halmd
