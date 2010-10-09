/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <cmath>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>


using namespace boost;
using namespace boost::fusion;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Assemble module options
 */
template <int dimension>
void core<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("backend",
#ifdef WITH_CUDA
         po::value<string>()->default_value("gpu"),
#else
         po::value<string>()->default_value("host"),
#endif
         "computing device type")
        ("dimension", po::value<int>()->default_value(3),
         "dimension of positional coordinates")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    using namespace lua_wrapper;
    register_any_converter<int>();
}

/**
 * Initialize simulation
 */
template <int dimension>
core<dimension>::core()
  // initialise attributes
  : step_counter_(0)
{
    LOG("dimension of positional coordinates: " << dimension);
    // FIXME LOG("MD simulation backend: " << vm["backend"].as<string>());
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void core<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * Prepare microscopic system state
 */
template <int dimension>
void core<dimension>::prepare()
{
    scoped_timer<timer> timer_(at_key<prepare_>(runtime_));
    particle->set();
    position->set();
    velocity->set();
    neighbour->update();
    force->compute();
}

/**
 * Perform a single MD integration step
 */
template <int dimension>
void core<dimension>::mdstep()
{
    scoped_timer<timer> timer_(at_key<mdstep_>(runtime_));
    integrator->integrate();
    if (neighbour->check()) {
        if (sort) {
            sort->order();
        }
        neighbour->update();
    }
    force->compute();
    integrator->finalize();

    // increment step counter
    step_counter_++;
}

template <typename T>
static void register_lua(lua_State* L, char const* class_name)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<T, shared_ptr<T> >(class_name)
                    .def(constructor<>())
                    .def("register_runtimes", &T::register_runtimes)
                    .def_readwrite("particle", &T::particle)
                    .def_readwrite("box", &T::box)
                    .def_readwrite("force", &T::force)
                    .def_readwrite("neighbour", &T::neighbour)
                    .def_readwrite("sort", &T::sort)
                    .def_readwrite("integrator", &T::integrator)
                    .def_readwrite("position", &T::position)
                    .def_readwrite("velocity", &T::velocity)
                    .property("step_counter", &T::step_counter)
                    .property("time", &T::time)
                    .def("prepare", &T::prepare)
                    .def("mdstep", &T::mdstep)
                    .scope
                    [
                        def("options", &T::options)
                    ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        bind(&register_lua<core<3> >, _1, "core_3_")
    ]
    [
        bind(&register_lua<core<2> >, _1, "core_2_")
    ];
}

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

} // namespace halmd
