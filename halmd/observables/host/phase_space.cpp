/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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
#include <functional>
#include <memory>

#include <halmd/observables/host/phase_space.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<particle_group_type> particle_group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<clock_type const> clock
  , std::shared_ptr<logger_type> logger
)
  : particle_(particle)
  , particle_group_(particle_group)
  , box_(box)
  , clock_(clock)
  , logger_(logger) {}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
std::shared_ptr<typename phase_space<dimension, float_type>::sample_type const>
phase_space<dimension, float_type>::acquire()
{
    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    group_array_type const& group = read_cache(particle_group_->ordered());
    position_array_type const& particle_position = read_cache(particle_->position());
    image_array_type const& particle_image = read_cache(particle_->image());
    velocity_array_type const& particle_velocity = read_cache(particle_->velocity());
    species_array_type const& particle_species = read_cache(particle_->species());
    mass_array_type const& particle_mass = read_cache(particle_->mass());

    LOG_TRACE("acquire sample");

    scoped_timer_type timer(runtime_.acquire);

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = std::make_shared<sample_type>(group.size(), clock_->step());
    }

    typename sample_type::position_array_type& sample_position = sample_->position();
    typename sample_type::velocity_array_type& sample_velocity = sample_->velocity();
    typename sample_type::species_array_type& sample_species = sample_->species();
    typename sample_type::mass_array_type& sample_mass = sample_->mass();

    // copy particle data using index map
    std::size_t tag = 0;
    for (std::size_t i : group) {
        // periodically extended particle position
        vector_type& r = sample_position[tag];
        r = particle_position[i];
        box_->extend_periodic(r, particle_image[i]);

        sample_velocity[tag] = particle_velocity[i];
        sample_species[tag] = particle_species[i];
        sample_mass[tag] = particle_mass[i];
        ++tag;
    }

    return sample_;
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::set(std::shared_ptr<sample_type const> sample)
{
    auto particle_position = make_cache_mutable(particle_->position());
    auto particle_image = make_cache_mutable(particle_->image());
    auto particle_velocity = make_cache_mutable(particle_->velocity());
    auto particle_species = make_cache_mutable(particle_->species());
    auto particle_mass = make_cache_mutable(particle_->mass());

    group_array_type const& group = read_cache(particle_group_->ordered());
    species_type nspecies = particle_->nspecies();

    typename sample_type::position_array_type const& sample_position = sample->position();
    typename sample_type::velocity_array_type const& sample_velocity = sample->velocity();
    typename sample_type::species_array_type const& sample_species = sample->species();
    typename sample_type::mass_array_type const& sample_mass = sample->mass();

    scoped_timer_type timer(runtime_.set);

    size_type tag = 0;
    for (size_type i : group) {
        vector_type& r = (*particle_position)[i];
        r = sample_position[tag];
        vector_type& image = (*particle_image)[i];
        image = 0;

        // The host implementation of reduce_periodic wraps the position at
        // most once around the box. This is more efficient during the
        // simulation. For setting arbitrary particle positions, however,
        // we must ensure that the final position is actually inside the
        // periodic box.
        vector_type shift;
        do {
            image += (shift = box_->reduce_periodic(r));
        } while (shift != vector_type(0));

        (*particle_velocity)[i] = sample_velocity[tag];
        species_type species = sample_species[tag];
        if (species >= nspecies) {
            throw std::invalid_argument("invalid species");
        }
        (*particle_species)[i] = species;
        (*particle_mass)[i] = sample_mass[tag];
        ++tag;
    }
}

template <typename phase_space_type>
static std::function<std::shared_ptr<typename phase_space_type::sample_type const> ()>
wrap_acquire(std::shared_ptr<phase_space_type> self)
{
    return [=]() {
        return self->acquire();
    };
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::position_array_type const& ()>
wrap_position(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::position_array_type const& {
        return self->acquire()->position();
    };
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::velocity_array_type const& ()>
wrap_velocity(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::velocity_array_type const& {
        return self->acquire()->velocity();
    };
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::species_array_type const& ()>
wrap_species(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::species_array_type const& {
        return self->acquire()->species();
    };
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::mass_array_type const& ()>
wrap_mass(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::mass_array_type const& {
        return self->acquire()->mass();
    };
}

template <typename phase_space_type>
static int wrap_dimension(phase_space_type const&)
{
    return phase_space_type::particle_type::vector_type::static_size;
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<phase_space>()
                    .property("acquire", &wrap_acquire<phase_space>)
                    .property("position", &wrap_position<phase_space>)
                    .property("velocity", &wrap_velocity<phase_space>)
                    .property("species", &wrap_species<phase_space>)
                    .property("mass", &wrap_mass<phase_space>)
                    .property("dimension", &wrap_dimension<phase_space>)
                    .def("set", &phase_space::set)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                            .def_readonly("reset", &runtime::reset)
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &phase_space::runtime_)

              , def("phase_space", &std::make_shared<phase_space
                   , std::shared_ptr<particle_type>
                   , std::shared_ptr<particle_group_type>
                   , std::shared_ptr<box_type const>
                   , std::shared_ptr<clock_type const>
                   , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    phase_space<3, double>::luaopen(L);
    phase_space<2, double>::luaopen(L);
#else
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#else
template class phase_space<3, float>;
template class phase_space<2, float>;
#endif

} // namespace observables
} // namespace host
} // namespace halmd
