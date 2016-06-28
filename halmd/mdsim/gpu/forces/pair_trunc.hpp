/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/forces/trunc/discontinuous.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

#include <memory>
#include <tuple>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type, typename trunc_type = mdsim::forces::trunc::discontinuous>
class pair_trunc
{
public:
    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;
    typedef neighbour neighbour_type;

    pair_trunc(
        std::shared_ptr<potential_type const> potential
      , std::shared_ptr<particle_type> particle1
      , std::shared_ptr<particle_type const> particle2
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<neighbour_type> neighbour
      , float_type aux_weight = 1
      , std::shared_ptr<trunc_type const> trunc = std::make_shared<trunc_type>()
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Check if the force cache (of the particle module) is up-to-date and if
     * not, mark the cache as dirty.
     */
    void check_cache();

    /**
     * Compute and apply the force to the particles in particle1.
     */
    void apply();

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::position_type position_type;
    typedef typename neighbour_type::array_type neighbour_array_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_trunc_wrapper<dimension, gpu_potential_type, trunc_type> gpu_wrapper;

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
    /** state of first system */
    std::shared_ptr<particle_type> particle1_;
    /** state of second system */
    std::shared_ptr<particle_type const> particle2_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** neighbour lists */
    std::shared_ptr<neighbour_type> neighbour_;
    /** weight for auxiliary variables */
    float_type aux_weight_;
    /** smoothing functor */
    std::shared_ptr<trunc_type const> trunc_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /** cache observer of force per particle */
    std::tuple<cache<>, cache<>> force_cache_;
    /** cache observer of auxiliary variables */
    std::tuple<cache<>, cache<>> aux_cache_;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
        accumulator_type compute_aux;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
pair_trunc<dimension, float_type, potential_type, trunc_type>::pair_trunc(
    std::shared_ptr<potential_type const> potential
  , std::shared_ptr<particle_type> particle1
  , std::shared_ptr<particle_type const> particle2
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<neighbour_type> neighbour
  , float_type aux_weight
  , std::shared_ptr<trunc_type const> trunc
  , std::shared_ptr<logger> logger
)
  : potential_(potential)
  , particle1_(particle1)
  , particle2_(particle2)
  , box_(box)
  , neighbour_(neighbour)
  , aux_weight_(aux_weight)
  , trunc_(trunc)
  , logger_(logger)
{
    if (std::min(potential_->size1(), potential_->size2()) < std::max(particle1_->nspecies(), particle2_->nspecies())) {
        throw std::invalid_argument("size of potential coefficients less than number of particle species");
    }
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
inline void pair_trunc<dimension, float_type, potential_type, trunc_type>::check_cache()
{
    cache<position_array_type> const& position1_cache = particle1_->position();
    cache<position_array_type> const& position2_cache = particle2_->position();

    auto current_state = std::tie(position1_cache, position2_cache);

    if (force_cache_ != current_state) {
        particle1_->mark_force_dirty();
    }

    if (aux_cache_ != current_state) {
        particle1_->mark_aux_dirty();
    }
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
inline void pair_trunc<dimension, float_type, potential_type, trunc_type>::apply()
{
    cache<position_array_type> const& position1_cache = particle1_->position();
    cache<position_array_type> const& position2_cache = particle2_->position();

    auto current_state = std::tie(position1_cache, position2_cache);

    if (particle1_->aux_enabled()) {
        compute_aux_();
        force_cache_ = current_state;
        aux_cache_ = force_cache_;
    }
    else {
        compute_();
        force_cache_ = current_state;
    }
    particle1_->force_zero_disable();
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
inline void pair_trunc<dimension, float_type, potential_type, trunc_type>::compute_()
{
    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    neighbour_array_type const& g_neighbour = read_cache(neighbour_->g_neighbour());
    auto force = make_cache_mutable(particle1_->mutable_force());

    LOG_TRACE("compute forces");

    scoped_timer_type timer(runtime_.compute);

    gpu_wrapper::kernel.r2.bind(position2);
    potential_->bind_textures();

    cuda::configure(particle1_->dim.grid, particle1_->dim.block);
    gpu_wrapper::kernel.compute(
        &*position1.begin()
      , &*force->begin()
      , &*g_neighbour.begin()
      , neighbour_->size()
      , neighbour_->stride()
      , nullptr
      , nullptr
      , particle1_->nspecies()
      , particle2_->nspecies()
      , static_cast<position_type>(box_->length())
      , *trunc_
      , particle1_->force_zero()
      , 1 // only relevant for kernel.compute_aux()
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
inline void pair_trunc<dimension, float_type, potential_type, trunc_type>::compute_aux_()
{
    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    neighbour_array_type const& g_neighbour = read_cache(neighbour_->g_neighbour());
    auto force = make_cache_mutable(particle1_->mutable_force());
    auto en_pot = make_cache_mutable(particle1_->mutable_potential_energy());
    auto stress_pot = make_cache_mutable(particle1_->mutable_stress_pot());

    LOG_TRACE("compute forces with auxiliary variables");

    scoped_timer_type timer(runtime_.compute_aux);

    gpu_wrapper::kernel.r2.bind(position2);
    potential_->bind_textures();

    float_type weight = aux_weight_;
    if (particle1_ == particle2_) {
        weight /= 2;
    }

    cuda::configure(particle1_->dim.grid, particle1_->dim.block);
    gpu_wrapper::kernel.compute_aux(
        &*position1.begin()
      , &*force->begin()
      , &*g_neighbour.begin()
      , neighbour_->size()
      , neighbour_->stride()
      , &*en_pot->begin()
      , &*stress_pot->begin()
      , particle1_->nspecies()
      , particle2_->nspecies()
      , static_cast<position_type>(box_->length())
      , *trunc_
      , particle1_->force_zero()
      , weight
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
void pair_trunc<dimension, float_type, potential_type, trunc_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                class_<pair_trunc>()
                    .def("check_cache", &pair_trunc::check_cache)
                    .def("apply", &pair_trunc::apply)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                            .def_readonly("compute_aux", &runtime::compute_aux)
                    ]
                    .def_readonly("runtime", &pair_trunc::runtime_)

              , def("pair_trunc", &std::make_shared<pair_trunc,
                    std::shared_ptr<potential_type const>
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<neighbour_type>
                  , float_type
                  , std::shared_ptr<trunc_type const>
                  , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP */
