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

#ifndef HALMD_MDSIM_HOST_FORCES_EXTERNAL_HPP
#define HALMD_MDSIM_HOST_FORCES_EXTERNAL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

#include <memory>
#include <tuple>

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

/**
 * template class for modules implementing short ranged potential forces
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
     * Compute and apply the force to the particles.
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
    typedef typename particle_type::species_array_type species_array_type;
    typedef typename particle_type::species_type species_type;
    typedef typename particle_type::size_type size_type;
    typedef typename particle_type::force_type force_type;
    typedef typename particle_type::stress_pot_type stress_pot_type;
    typedef typename particle_type::en_pot_type en_pot_type;

    /** compute forces */
    void compute_();
    /** compute forces with auxiliary variables */
    void compute_aux_();

    /** pair potential */
    std::shared_ptr<potential_type const> potential_;
    /** state of particle system */
    std::shared_ptr<particle_type> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /** cache observer of force: (position, species) */
    std::tuple<cache<>, cache<>> force_cache_;
    /** cache observer of auxiliary variables: (position, species) */
    std::tuple<cache<>, cache<>> aux_cache_;

    /** store signal connections */
    signal_type on_prepend_apply_;
    signal_type on_append_apply_;

    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        utility::profiler::accumulator_type compute;
        utility::profiler::accumulator_type compute_aux;
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
    cache<species_array_type> const& species_cache = particle_->species();

    auto current_state = std::tie(position_cache, species_cache);

    if (force_cache_ != current_state) {
        particle_->mark_force_dirty();
    }

    if (aux_cache_ != current_state) {
        particle_->mark_aux_dirty();
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void external<dimension, float_type, potential_type>::apply()
{
    // process slot functions associated with signal
    on_prepend_apply_();

    cache<position_array_type> const& position_cache = particle_->position();
    cache<species_array_type> const& species_cache = particle_->species();

    auto current_state = std::tie(position_cache, species_cache);

    if (particle_->aux_enabled()) {
        compute_aux_();
        force_cache_ = current_state;
        aux_cache_ = force_cache_;
    }
    else {
        compute_();
        force_cache_ = current_state;
    }
    particle_->force_zero_disable();

    // process slot functions associated with signal
    on_append_apply_();
}

template <int dimension, typename float_type, typename potential_type>
inline void external<dimension, float_type, potential_type>::compute_()
{
    auto force = make_cache_mutable(particle_->mutable_force());

    position_array_type const& position = read_cache(particle_->position());
    species_array_type const& species   = *particle_->species();
    size_type nparticle = particle_->nparticle();

    LOG_DEBUG("compute forces");

    scoped_timer_type timer(runtime_.compute);

    // reset the force and auxiliary variables to zero if necessary
    if (particle_->force_zero()) {
        std::fill(force->begin(), force->end(), 0);
    }

    for (size_type i = 0; i < nparticle; ++i) {
        // reduced particle position
        position_type r = position[i];
        box_->reduce_periodic(r);
        // particle species
        species_type s = species[i];

        // evaluate potential
        force_type f;
        en_pot_type en_pot_;
        std::tie(f, en_pot_) = (*potential_)(r, s);

        // add force contribution
        (*force)[i] += f;
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void external<dimension, float_type, potential_type>::compute_aux_()
{
    auto force      = make_cache_mutable(particle_->mutable_force());
    auto en_pot     = make_cache_mutable(particle_->mutable_potential_energy());
    auto stress_pot = make_cache_mutable(particle_->mutable_stress_pot());

    position_array_type const& position = read_cache(particle_->position());
    species_array_type const& species   = *particle_->species();
    size_type nparticle = particle_->nparticle();

    LOG_DEBUG("compute forces with auxiliary variables");

    scoped_timer_type timer(runtime_.compute_aux);

    // reset the force and auxiliary variables to zero if necessary
    if (particle_->force_zero()) {
        std::fill(force->begin(), force->end(), 0);
        std::fill(en_pot->begin(), en_pot->end(), 0);
        std::fill(stress_pot->begin(), stress_pot->end(), 0);
    }

    for (size_type i = 0; i < nparticle; ++i) {
        // reduced particle position
        position_type r = position[i];
        box_->reduce_periodic(r);
        // particle species
        species_type s = species[i];

        // evaluate potential
        force_type f;
        en_pot_type en_pot_;
        std::tie(f, en_pot_) = (*potential_)(r, s);

        // add force contribution
        (*force)[i] += f;

        // add contribution to potential energy
        (*en_pot)[i] += en_pot_;

        // add potential part of stress tensor
//        (*stress_pot)[i] += make_stress_tensor(r, f); FIXME
    }
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
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_EXTERNAL_HPP */
