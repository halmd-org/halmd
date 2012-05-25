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
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>

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
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<particle_group_type const> particle_group
  , boost::shared_ptr<box_type const> box
  , boost::shared_ptr<clock_type const> clock
  , boost::shared_ptr<logger_type> logger
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
boost::shared_ptr<typename phase_space<dimension, float_type>::sample_type const>
phase_space<dimension, float_type>::acquire()
{
    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    scoped_timer_type timer(runtime_.acquire);

    LOG_TRACE("acquire sample");

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = boost::make_shared<sample_type>(particle_group_->size(), clock_->step());
    }

    typename particle_type::position_array_type& particle_position = particle_->position();
    typename particle_type::image_array_type& particle_image = particle_->image();
    typename particle_type::velocity_array_type& particle_velocity = particle_->velocity();
    typename particle_type::species_array_type& particle_species = particle_->species();

    typename sample_type::position_array_type& sample_position = sample_->position();
    typename sample_type::velocity_array_type& sample_velocity = sample_->velocity();
    typename sample_type::species_array_type& sample_species = sample_->species();

    // copy particle data using index map
    std::size_t tag = 0;
    for (std::size_t i : *particle_group_) {
        assert(i < particle_->nparticle());

        // periodically extended particle position
        vector_type& r = sample_position[tag];
        r = particle_position[i];
        box_->extend_periodic(r, particle_image[i]);

        sample_velocity[tag] = particle_velocity[i];
        sample_species[tag] = particle_species[i];
        ++tag;
    }

    return sample_;
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::set(boost::shared_ptr<sample_type const> sample)
{
    scoped_timer_type timer(runtime_.set);

    typename particle_type::position_array_type& particle_position = particle_->position();
    typename particle_type::image_array_type& particle_image = particle_->image();
    typename particle_type::velocity_array_type& particle_velocity = particle_->velocity();
    typename particle_type::species_array_type& particle_species = particle_->species();

    typename sample_type::position_array_type const& sample_position = sample->position();
    typename sample_type::velocity_array_type const& sample_velocity = sample->velocity();
    typename sample_type::species_array_type const& sample_species = sample->species();

    std::size_t tag = 0;
    for (std::size_t i : *particle_group_) {
        assert(i < particle_->nparticle());
        vector_type& r = particle_position[i];
        r = sample_position[tag];
        vector_type& image = particle_image[i];
        image = 0;

        // The host implementation of reduce_periodic wraps the position at
        // most once around the box. This is more efficient during the
        // simulation. For setting arbitrary particle positions, however,
        // we must ensure that the final position is actually inside the
        // periodic box.
        vector_type shift;
        do {
            image += (shift = box_->reduce_periodic(r));
        } while(shift != vector_type(0));

        particle_velocity[i] = sample_velocity[tag];
        particle_species[i] = sample_species[tag];
        ++tag;
    }
}

template <typename phase_space_type, typename sample_type>
static boost::function<boost::shared_ptr<sample_type const> ()>
wrap_acquire(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&phase_space_type::acquire, phase_space);
}

template <typename phase_space_type>
static typename phase_space_type::sample_type::position_array_type const&
position(boost::shared_ptr<phase_space_type> const& phase_space)
{
    return phase_space->acquire()->position();
}

template <typename phase_space_type>
static boost::function<typename phase_space_type::sample_type::position_array_type const& ()>
wrap_position(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&position<phase_space_type>, phase_space);
}

template <typename phase_space_type>
static typename phase_space_type::sample_type::velocity_array_type const&
velocity(boost::shared_ptr<phase_space_type> const& phase_space)
{
    return phase_space->acquire()->velocity();
}

template <typename phase_space_type>
static boost::function<typename phase_space_type::sample_type::velocity_array_type const& ()>
wrap_velocity(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&velocity<phase_space_type>, phase_space);
}

/**
 * Translate species from internal (0-based) to external (1-based) representation.
 *
 * This function returns a const reference to a species array that is
 * stored in the functor, and passed by reference to this function.
 */
template <typename phase_space_type>
static typename phase_space_type::sample_type::species_array_type const&
species(boost::shared_ptr<phase_space_type> const& phase_space, boost::shared_ptr<typename phase_space_type::sample_type::species_array_type>& species)
{
    boost::shared_ptr<typename phase_space_type::sample_type const> sample = phase_space->acquire();
    species.reset(new typename phase_space_type::sample_type::species_array_type(sample->species().size()));
    std::transform(
        sample->species().begin()
      , sample->species().end()
      , species->begin()
      , [](typename phase_space_type::sample_type::species_array_type::value_type s) {
            return s + 1;
        }
    );
    return *species;
}

template <typename phase_space_type>
static boost::function<typename phase_space_type::sample_type::species_array_type const& ()>
wrap_species(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&species<phase_space_type>, phase_space, boost::shared_ptr<typename phase_space_type::sample_type::species_array_type>());
}

template <typename phase_space_type>
static int wrap_dimension(phase_space_type const&)
{
    return phase_space_type::particle_type::vector_type::static_size;
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<phase_space>()
                    .property("acquire", &wrap_acquire<phase_space, sample_type>)
                    .property("position", &wrap_position<phase_space>)
                    .property("velocity", &wrap_velocity<phase_space>)
                    .property("species", &wrap_species<phase_space>)
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

              , def("phase_space", &boost::make_shared<phase_space
                   , boost::shared_ptr<particle_type>
                   , boost::shared_ptr<particle_group_type const>
                   , boost::shared_ptr<box_type const>
                   , boost::shared_ptr<clock_type const>
                   , boost::shared_ptr<logger_type>
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
