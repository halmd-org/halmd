/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/forces/pair_full_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

#include <memory>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_full
{
public:
    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;
    typedef logger logger_type;

    pair_full(
        std::shared_ptr<potential_type const> potential
      , std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /**
     * Check if the force cache (of the particle module) is up-to-date and if
     * not, mark the cache as dirty.
     */
    void check_cache();

    /**
     * Compute and apply the force to the particles.
     */
    void apply();


    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::position_type position_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_full_wrapper<dimension, gpu_potential_type> gpu_wrapper;

    typedef typename particle_type::force_array_type force_array_type;
    typedef typename particle_type::en_pot_array_type en_pot_array_type;
    typedef typename particle_type::stress_pot_type stress_pot_type;
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;

    /** compute forces */
    void compute_();
    /** compute forces with auxiliary variables */
    void compute_aux_();

    /** pair potential */
    std::shared_ptr<potential_type const> potential_;
    /** system state */
    std::shared_ptr<particle_type> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** cache observer of force per particle */
    cache<> force_cache_;
    /** cache observer of auxiliary variables */
    cache<> aux_cache_;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type>
pair_full<dimension, float_type, potential_type>::pair_full(
    std::shared_ptr<potential_type const> potential
  , std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger_type> logger
)
  : potential_(potential)
  , particle_(particle)
  , box_(box)
  , logger_(logger)
{
    if (std::min(potential_->size1(), potential_->size2()) < particle_->nspecies()) {
        throw std::invalid_argument("size of potential coefficients less than number of particle species");
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::check_cache()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (force_cache_ != position_cache ||
        (particle_->aux_enabled() && aux_cache_ != position_cache)) {
        particle_->mark_force_dirty();
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::apply()
{
    if (particle_->aux_enabled()) {
        compute_aux_();
        force_cache_ = particle_->position();
        aux_cache_ = force_cache_;
    }
    else {
        compute_();
        force_cache_ = particle_->position();
    }
    particle_->force_zero_disable();
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::compute_()
{
    LOG_TRACE("compute forces");

    position_array_type const& position = read_cache(particle_->position());
    auto force = make_cache_mutable(particle_->mutable_force());

    scoped_timer_type timer(runtime_.compute);

    potential_->bind_textures();

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    gpu_wrapper::kernel.compute(
        &*force->begin()
        , &*position.begin()
        , nullptr
        , nullptr
        , particle_->nparticle()
        , particle_->nspecies()
        , particle_->nspecies()
        , static_cast<position_type>(box_->length())
        , particle_->force_zero()
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::compute_aux_()
{
    LOG_TRACE("compute forces with auxiliary variables");

    position_array_type const& position = read_cache(particle_->position());
    auto force = make_cache_mutable(particle_->mutable_force());
    auto en_pot = make_cache_mutable(particle_->mutable_potential_energy());
    auto stress_pot = make_cache_mutable(particle_->mutable_stress_pot());

    scoped_timer_type timer(runtime_.compute);

    potential_->bind_textures();

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    gpu_wrapper::kernel.compute_aux(
        &*force->begin()
        , &*position.begin()
        , &*en_pot->begin()
        , &*stress_pot->begin()
        , particle_->nparticle()
        , particle_->nspecies()
        , particle_->nspecies()
        , static_cast<position_type>(box_->length())
        , particle_->force_zero()
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
void pair_full<dimension, float_type, potential_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                class_<pair_full>()
                    .def("check_cache", &pair_full::check_cache)
                    .def("apply", &pair_full::apply)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &pair_full::runtime_)

              , def("pair_full", &std::make_shared<pair_full,
                    std::shared_ptr<potential_type const>
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP */
