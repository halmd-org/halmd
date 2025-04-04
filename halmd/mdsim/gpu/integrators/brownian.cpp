/*
 * Copyright © 2011-2024 Felix Höfling
 * Copyright © 2011      Michael Kopp
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

#include <boost/numeric/ublas/io.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/integrators/brownian.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type, typename RandomNumberGenerator>
brownian<dimension, float_type, RandomNumberGenerator>::brownian(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , std::shared_ptr<box_type const> box
  , double timestep
  , double temperature
  , scalar_container_type const& diffusion
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , random_(random)
  , box_(box)
  , diffusion_(diffusion)
  , g_param_(particle->nspecies())
  , logger_(logger)
{
    if (diffusion_.size() != particle_->nspecies()) {
        throw std::invalid_argument("diffusion constants have mismatching shape");
    }

    set_timestep(timestep);
    set_temperature(temperature);       // assigns g_param_: noise, mobility

    LOG("diffusion constants: " << diffusion_);
}

/**
 * set integration timestep
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::set_timestep(double timestep)
{
    timestep_ = timestep;
    LOG("integration timestep: " << float(timestep_));
}

/**
 * set temperature of the heat bath
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::set_temperature(double temperature)
{
    temperature_ = temperature;
    LOG("temperature: " << float(temperature_));

    // re-compute and copy parameters to CUDA device
    cuda::memory::host::vector<float2> param(g_param_.size());
    scalar_container_type noise(param.size());
    scalar_container_type mobility(param.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 2> p(0);
        p[brownian_param::NOISE]    = noise[i] = sqrt(2 * diffusion_[i]);
        p[brownian_param::MOBILITY] = mobility[i] = diffusion_[i] / temperature_;
        param[i] = p;
    }
    cuda::copy(param.begin(), param.end(), g_param_.begin());

    LOG_INFO("noise strengths: " << noise);
    LOG_INFO("mobility constants: " << mobility);
}

/**
 * perform Brownian integration: update positions from random distribution
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::integrate()
{
    LOG_TRACE("update positions")

    force_array_type const& force = read_cache(particle_->force());

    // invalidate the particle caches only after accessing the force!
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());

    scoped_timer_type timer(runtime_.integrate);

    try {
        cuda::texture<float2> t_param(g_param_);

        // use CUDA execution dimensions of 'random' since
        // the kernel makes use of the random number generator
        wrapper_type::kernel.integrate.configure(random_->rng().dim.grid,
            random_->rng().dim.block);
        wrapper_type::kernel.integrate(
            t_param
          , position->data()
          , image->data()
          , force.data()
          , timestep_
          , particle_->nparticle()
          , static_cast<vector_type>(box_->length())
          , random_->rng().rng()
        );
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to stream Brownian integration on GPU");
        throw;
    }
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void brownian<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<brownian>()
                    .def("integrate", &brownian::integrate)
                    .def("set_timestep", &brownian::set_timestep)
                    .def("set_temperature", &brownian::set_temperature)
                    .property("timestep", &brownian::timestep)
                    .property("temperature", &brownian::temperature)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                    ]
                    .def_readonly("runtime", &brownian::runtime_)

              , def("brownian", &std::make_shared<brownian
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<random_type>
                  , std::shared_ptr<box_type const>
                  , double
                  , double
                  , scalar_container_type const&
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_integrators_brownian(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    brownian<2, float, halmd::random::gpu::rand48>::luaopen(L);
    brownian<3, float, halmd::random::gpu::rand48>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    brownian<2, dsfloat, halmd::random::gpu::rand48>::luaopen(L);
    brownian<3, dsfloat, halmd::random::gpu::rand48>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class brownian<2, float, halmd::random::gpu::rand48>;
template class brownian<3, float, halmd::random::gpu::rand48>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class brownian<2, dsfloat, halmd::random::gpu::rand48>;
template class brownian<3, dsfloat, halmd::random::gpu::rand48>;
#endif

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
