/*
 * Copyright © 2010-2014 Felix Höfling
 * Copyright © 2016      Sutapa Roy
 * Copyright © 2013-2014 Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_FORCES_EXTERNAL_HPP
#define HALMD_MDSIM_GPU_FORCES_EXTERNAL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/forces/external_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

#include <memory>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing forces due to an external potential
 */
template <int dimension, typename float_type, typename potential_type>
class external
{
public:
    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;

    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    external(
        std::shared_ptr<potential_type const> potential
      , std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
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
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef external_wrapper<dimension, gpu_potential_type> gpu_wrapper;

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
    std::shared_ptr<particle_type> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /** cache observer of force on particle */
    cache<> force_cache_;
    /** cache observer of auxiliary variables */
    cache<> aux_cache_;

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
external<dimension, float_type, potential_type>::external(
    std::shared_ptr<potential_type const> potential
  , std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger> logger
)
  : potential_(potential)
  , particle_(particle)
  , box_(box)
  , logger_(logger)
{
    if (potential_->size() < particle_->nspecies()) {
        throw std::invalid_argument("size of potential coefficients less than number of particle species");
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void external<dimension, float_type, potential_type>::check_cache()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (force_cache_ != position_cache) {
        particle_->mark_force_dirty();
    }

    if (aux_cache_ != position_cache) {
        particle_->mark_aux_dirty();
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void external<dimension, float_type, potential_type>::apply()
{
    // process slot functions associated with signal
    on_prepend_apply_();

    cache<position_array_type> const& position_cache = particle_->position();

    if (particle_->aux_enabled()) {
        compute_aux_();
        force_cache_ = position_cache;
        aux_cache_ = force_cache_;
    }
    else {
        compute_();
        force_cache_ = position_cache;
    }
    particle_->force_zero_disable();

    // process slot functions associated with signal
    on_append_apply_();
}

template <int dimension, typename float_type, typename potential_type>
inline void external<dimension, float_type, potential_type>::compute_()
{
    position_array_type const& position = read_cache(particle_->position());
    auto force = make_cache_mutable(particle_->mutable_force());

    LOG_DEBUG("compute forces");

    scoped_timer_type timer(runtime_.compute);

    configure_kernel(gpu_wrapper::kernel.compute, particle_->dim(), true);
    gpu_wrapper::kernel.compute(
        potential_->get_gpu_potential()
      , position.data()
      , &*force->begin()
      , nullptr
      , nullptr
      , static_cast<position_type>(box_->length())
      , particle_->force_zero()
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
inline void external<dimension, float_type, potential_type>::compute_aux_()
{
    position_array_type const& position = read_cache(particle_->position());
    auto force = make_cache_mutable(particle_->mutable_force());
    auto en_pot = make_cache_mutable(particle_->mutable_potential_energy());
    auto stress_pot = make_cache_mutable(particle_->mutable_stress_pot());

    LOG_DEBUG("compute forces with auxiliary variables");

    scoped_timer_type timer(runtime_.compute_aux);

    configure_kernel(gpu_wrapper::kernel.compute_aux, particle_->dim(), true);
    gpu_wrapper::kernel.compute_aux(
        potential_->get_gpu_potential()
      , position.data()
      , &*force->begin()
      , &*en_pot->begin()
      , &*stress_pot->begin()
      , static_cast<position_type>(box_->length())
      , particle_->force_zero()
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
void external<dimension, float_type, potential_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                class_<external>()
                    .def("check_cache", &external::check_cache)
                    .def("apply", &external::apply)
                    .def("on_prepend_apply", &external::on_prepend_apply)
                    .def("on_append_apply", &external::on_append_apply)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                            .def_readonly("compute_aux", &runtime::compute_aux)
                    ]
                    .def_readonly("runtime", &external::runtime_)

              , def("external", &std::make_shared<external,
                    std::shared_ptr<potential_type const>
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
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

#endif /* ! HALMD_MDSIM_GPU_FORCES_EXTERNAL_HPP */
