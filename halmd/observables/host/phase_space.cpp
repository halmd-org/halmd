/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <boost/make_shared.hpp>

#include <halmd/observables/host/phase_space.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    boost::shared_ptr<sample_type> sample
  , boost::shared_ptr<particle_type const> particle
  , boost::shared_ptr<box_type const> box
  , boost::shared_ptr<clock_type const> clock
  , boost::shared_ptr<logger_type> logger
)
  : sample_(sample)
  , particle_(particle)
  , box_(box)
  , clock_(clock)
  , logger_(logger)
{
}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
void phase_space<dimension, float_type>::acquire()
{
    scoped_timer_type timer(runtime_.acquire);

    if (sample_->step == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return;
    }

    LOG_TRACE("acquire sample");

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_->reset();
    }

    for (size_t i = 0; i < particle_->nbox; ++i) {
        unsigned int type = particle_->type[i];
        unsigned int tag = particle_->tag[i];

        // periodically extended particle position
        assert(type < sample_->r.size());
        assert(tag < sample_->r[type]->size());
        vector_type& r = (*sample_->r[type])[tag] = particle_->r[i];
        box_->extend_periodic(r, particle_->image[i]);

        // particle velocity
        assert(type < sample_->v.size());
        assert(tag < sample_->v[type]->size());
        (*sample_->v[type])[tag] = particle_->v[i];
    }
    sample_->step = clock_->step();
}

template <int dimension, typename float_type>
static int wrap_dimension(phase_space<dimension, float_type> const&)
{
    return dimension;
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
            class_<phase_space, boost::shared_ptr<_Base>, _Base>(class_name.c_str())
                .property("dimension", &wrap_dimension<dimension, float_type>)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                        .def_readonly("reset", &runtime::reset)
                ]
                .def_readonly("runtime", &phase_space::runtime_)

          , def("phase_space", &boost::make_shared<phase_space
               , boost::shared_ptr<sample_type>
               , boost::shared_ptr<particle_type const>
               , boost::shared_ptr<box_type const>
               , boost::shared_ptr<clock_type const>
               , boost::shared_ptr<logger_type>
            >)
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
