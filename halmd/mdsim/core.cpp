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
        ("dimension", po::value<int>()->default_value(3),
         "dimension of positional coordinates")
        ;
}

template <int dimension>
void core<dimension>::write_parameters(H5::Group const& param) const
{
    h5xx::write_attribute(param, "dimension", dimension);
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
    LOG_TRACE("performing MD step #" << step_counter_);
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

/**
 * wrap dimension template parameter for Lua
 */
template <int dimension>
static int get_dimension(core<dimension> const&)
{
    return dimension;
}

template <int dimension>
void core<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    string class_name("core_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<core, shared_ptr<core> >(class_name.c_str())
                    .def(constructor<>())
                    .def("register_runtimes", &core::register_runtimes)
                    .def("write_parameters", &core::write_parameters)
                    .def_readwrite("particle", &core::particle)
                    .def_readwrite("box", &core::box)
                    .def_readwrite("force", &core::force)
                    .def_readwrite("neighbour", &core::neighbour)
                    .def_readwrite("sort", &core::sort)
                    .def_readwrite("integrator", &core::integrator)
                    .def_readwrite("position", &core::position)
                    .def_readwrite("velocity", &core::velocity)
                    .property("dimension", &get_dimension<dimension>)
                    .property("step_counter", &core::step_counter)
                    .def("prepare", &core::prepare)
                    .def("mdstep", &core::mdstep)
                    .scope
                    [
                        def("options", &core::options)
                    ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &core<3>::luaopen
    ]
    [
        &core<2>::luaopen
    ];
}

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

} // namespace halmd
