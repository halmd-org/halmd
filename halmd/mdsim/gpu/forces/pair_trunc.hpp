/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2016      Sutapa Roy
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2012 Peter Colberg
 * Copyright © 2021      Jaslo Ziska
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
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>
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
template <int dimension, typename float_type, typename potential_type>
class pair_trunc
{
public:
    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;
    typedef neighbour neighbour_type;
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    pair_trunc(
        std::shared_ptr<potential_type> potential
      , std::shared_ptr<particle_type> particle1
      , std::shared_ptr<particle_type const> particle2
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<neighbour_type> neighbour
      , float aux_weight = 1
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
     * Connect slot functions to signals
     */
    connection on_prepend_apply(slot_function_type const& slot)
    {
        return on_prepend_apply_.connect(slot);
    }

    connection on_append_apply(slot_function_type const& slot)
    {
        return on_append_apply_.connect(slot);
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::position_type position_type;
    typedef typename neighbour_type::array_type neighbour_array_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_trunc_wrapper<dimension, gpu_potential_type> gpu_wrapper;

    typedef typename particle_type::force_array_type force_array_type;
    typedef typename particle_type::en_pot_array_type en_pot_array_type;
    typedef typename particle_type::stress_pot_type stress_pot_type;
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;

    /** compute forces */
    void compute_();
    /** compute forces with auxiliary variables */
    void compute_aux_();

    /** pair potential */
    std::shared_ptr<potential_type> potential_;
    /** state of first system */
    std::shared_ptr<particle_type> particle1_;
    /** state of second system */
    std::shared_ptr<particle_type const> particle2_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** neighbour lists */
    std::shared_ptr<neighbour_type> neighbour_;
    /** weight for auxiliary variables */
    float aux_weight_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /** cache observer of force per particle */
    std::tuple<cache<>, cache<>> force_cache_;
    /** cache observer of auxiliary variables */
    std::tuple<cache<>, cache<>> aux_cache_;

    /** store signal connections */
    signal_type on_prepend_apply_;
    signal_type on_append_apply_;

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

template <int dimension, typename float_type, typename potential_type>
pair_trunc<dimension, float_type, potential_type>::pair_trunc(
    std::shared_ptr<potential_type> potential
  , std::shared_ptr<particle_type> particle1
  , std::shared_ptr<particle_type const> particle2
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<neighbour_type> neighbour
  , float aux_weight
  , std::shared_ptr<logger> logger
)
  : potential_(potential)
  , particle1_(particle1)
  , particle2_(particle2)
  , box_(box)
  , neighbour_(neighbour)
  , aux_weight_(aux_weight)
  , logger_(logger)
{
    if (std::min(potential_->size1(), potential_->size2()) < std::max(particle1_->nspecies(), particle2_->nspecies())) {
        throw std::invalid_argument("size of potential coefficients less than number of particle species");
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_trunc<dimension, float_type, potential_type>::check_cache()
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

template <int dimension, typename float_type, typename potential_type>
inline void pair_trunc<dimension, float_type, potential_type>::apply()
{
    // process slot functions associated with signal
    on_prepend_apply_();

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
    // process slot functions associated with signal
    on_append_apply_();

}

template <int dimension, typename float_type, typename potential_type>
inline void pair_trunc<dimension, float_type, potential_type>::compute_()
{
    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    neighbour_array_type const& g_neighbour = read_cache(neighbour_->g_neighbour());
    auto force = make_cache_mutable(particle1_->mutable_force());

    cuda::texture<float4> t_r2(position2);

    LOG_TRACE("compute forces");

    scoped_timer_type timer(runtime_.compute);

    configure_kernel(
        gpu_wrapper::kernel.compute
      , particle1_->array_size() * gpu_wrapper::kernel.nparallel_particles
      , true
    );
    gpu_wrapper::kernel.compute(
        potential_->get_gpu_potential()
      , position1.data()
      , t_r2
      , force->data()
      , g_neighbour.data()
      , neighbour_->size()
      , nullptr
      , nullptr
      , particle1_->nspecies()
      , particle2_->nspecies()
      , static_cast<position_type>(box_->length())
      , particle1_->force_zero()
      , 1 // only relevant for kernel.compute_aux()
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_trunc<dimension, float_type, potential_type>::compute_aux_()
{
    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    neighbour_array_type const& g_neighbour = read_cache(neighbour_->g_neighbour());
    auto force = make_cache_mutable(particle1_->mutable_force());
    auto en_pot = make_cache_mutable(particle1_->mutable_potential_energy());
    auto stress_pot = make_cache_mutable(particle1_->mutable_stress_pot());

    cuda::texture<float4> t_r2(position2);

    LOG_TRACE("compute forces with auxiliary variables");

    scoped_timer_type timer(runtime_.compute_aux);

    float weight = aux_weight_;
    if (particle1_ == particle2_) {
        weight /= 2;
    }

    configure_kernel(
        gpu_wrapper::kernel.compute_aux
      , particle1_->array_size() * gpu_wrapper::kernel.nparallel_particles
      , true
    );
    gpu_wrapper::kernel.compute_aux(
        potential_->get_gpu_potential()
      , position1.data()
      , t_r2
      , force->data()
      , g_neighbour.data()
      , neighbour_->size()
      , &*en_pot->begin()
      , &*stress_pot->begin()
      , particle1_->nspecies()
      , particle2_->nspecies()
      , static_cast<position_type>(box_->length())
      , particle1_->force_zero()
      , weight
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::luaopen(lua_State* L)
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
                    .def("on_prepend_apply", &pair_trunc::on_prepend_apply)
                    .def("on_append_apply", &pair_trunc::on_append_apply)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                            .def_readonly("compute_aux", &runtime::compute_aux)
                    ]
                    .def_readonly("runtime", &pair_trunc::runtime_)

              , def("pair_trunc", &std::make_shared<pair_trunc,
                    std::shared_ptr<potential_type>
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<neighbour_type>
                  , float
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
