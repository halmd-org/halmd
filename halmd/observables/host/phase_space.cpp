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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <halmd/observables/host/phase_space.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    shared_ptr<particle_group_type const> particle_group
  , shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , shared_ptr<logger_type> logger
)
  : particle_group_(particle_group)
  , particle_(particle_group->particle())
  , box_(box)
  , clock_(clock)
  , logger_(logger)
{
}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
shared_ptr<typename phase_space<dimension, float_type>::sample_type const>
phase_space<dimension, float_type>::acquire()
{
    scoped_timer_type timer(runtime_.acquire);

    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    LOG_TRACE("acquire sample");

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = make_shared<sample_type>(particle_group_->size(), clock_->step());
    }

    typename particle_type::position_array_type& particle_position = particle_->r;
    typename particle_type::image_array_type& particle_image = particle_->image;
    typename particle_type::velocity_array_type& particle_velocity = particle_->v;
    typename particle_type::species_array_type& particle_species = particle_->type;

    typename sample_type::position_array_type& sample_position = sample_->position();
    typename sample_type::velocity_array_type& sample_velocity = sample_->velocity();
    typename sample_type::species_array_type& sample_species = sample_->species();

    // copy particle data using index map
    typename particle_group_type::map_iterator idx = particle_group_->map();
    for (unsigned int i = 0; i < particle_group_->size(); ++i, ++idx) {

        assert(*idx < particle_->nbox);

        // periodically extended particle position
        vector_type& r = sample_position[i] = particle_position[*idx];
        box_->extend_periodic(r, particle_image[*idx]);

        sample_velocity[i] = particle_velocity[*idx];
        sample_species[i] = particle_species[*idx];
    }

    return sample_;
}

template <typename phase_space_type, typename sample_type>
static function<shared_ptr<sample_type const> ()>
wrap_acquire(shared_ptr<phase_space_type> phase_space)
{
    return bind(&phase_space_type::acquire, phase_space);
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
    static string const class_name(demangled_name<phase_space>());
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space>(class_name.c_str())
                .property("acquire", &wrap_acquire<phase_space, sample_type>)
                .property("dimension", &wrap_dimension<phase_space>)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                        .def_readonly("reset", &runtime::reset)
                ]
                .def_readonly("runtime", &phase_space::runtime_)

          , namespace_("host")
            [
                def("phase_space", &make_shared<phase_space
                   , shared_ptr<particle_group_type const>
                   , shared_ptr<box_type const>
                   , shared_ptr<clock_type const>
                   , shared_ptr<logger_type>
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
