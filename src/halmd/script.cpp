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

#include <algorithm>
#include <boost/bind.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/lua.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd
{

/**
 * Assemble module options
 */
template <int dimension>
void script<dimension>::options(po::options_description& desc)
{
    po::options_description group("Simulation");
    group.add_options()
        ("steps,s", po::value<uint64_t>()->default_value(10000),
         "number of simulation steps")
        ("time,t", po::value<double>(),
         "total simulation time")
        ;
    desc.add(group);
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void script<dimension>::depends()
{
    modules::depends<_Self, core_type>::required();
    modules::depends<_Self, sampler_type>::required();
    modules::depends<_Self, profile_writer_type>::required();
}

template <int dimension>
script<dimension>::script(modules::factory& factory, po::options const& vm)
  : _Base(factory, vm)
  // dependency injection
  , core(modules::fetch<core_type>(factory, vm))
  , sampler(modules::fetch<sampler_type>(factory, vm))
  , profile_writers(modules::fetch<profile_writer_type>(factory, vm))
  , L_(luaL_newstate(), lua_close) //< create Lua state
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    luaL_openlibs(L); //< load Lua standard libraries

    // parse options
    if (vm["steps"].defaulted() && !vm["time"].empty()) {
        time_ = vm["time"].as<double>();
        steps_ = static_cast<uint64_t>(round(time_ / core->integrator->timestep()));
    }
    else {
        steps_ = vm["steps"].as<uint64_t>();
        time_ = steps_ * core->integrator->timestep();
    }

    LOG("number of integration steps: " << steps_);
    LOG("integration time: " << time_);
}

/**
 * Load HALMD Lua wrapper
 *
 * Register C++ classes with Lua.
 */
template <int dimension>
void script<dimension>::load_wrapper()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;
    using luabind::module; //< FIXME namespace conflicts

    open(L); //< setup global structures and Lua class support

    for_each(
        lua_registry::get()->begin()
      , lua_registry::get()->end()
      , bind(&module_::operator[], module(L), _1)
    );
}

/**
 * Load HALMD Lua library
 */
template <int dimension>
void script<dimension>::load_library()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    string path;
    path.append( HALMD_BINARY_DIR "/lib/?.lua" ";" );
    path.append( HALMD_BINARY_DIR "/lib/?/init.lua" ";" );
    path.append( HALMD_SOURCE_DIR "/lib/?.lua" ";" );
    path.append( HALMD_SOURCE_DIR "/lib/?/init.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/share/?.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/share/?/init.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/lib/?.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/lib/?/init.lua" ";" );
    path.append( object_cast<string>(globals(L)["package"]["path"]) );
    globals(L)["package"]["path"] = path;

    int status = luaL_dostring(L, "require('halmd')");

    if (status != 0) {
        LOG_ERROR("[Lua] " << lua_tostring(L, -1));
        lua_pop(L, 1); //< remove error message
        throw std::runtime_error("Lua error");
    }
}

/**
 * Run simulation
 */
template <int dimension>
void script<dimension>::run()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    LOG("starting simulation");

    int status = luaL_dostring(L, "halmd.lib.run()");

    if (status != 0) {
        LOG_ERROR("[Lua] " << lua_tostring(L, -1));
        lua_pop(L, 1); //< remove error message
        throw std::runtime_error("Lua error");
    }

    LOG("finished simulation");

    core->prepare();
    sampler->sample(true);

    LOG("starting NVE ensemble run");

    while (core->step_counter() < steps_) {
        // perform complete MD integration step
        core->mdstep();

        // sample system state and properties,
        // force sampling after last integration step
        sampler->sample(core->step_counter() == steps_);
    }

    LOG("finished NVE ensemble run");

    for_each(
        profile_writers.begin()
      , profile_writers.end()
      , bind(&profile_writer_type::write, _1)
    );
}

// explicit instantiation
template class script<3>;
template class script<2>;

template class module<script<3> >;
template class module<script<2> >;

} // namespace halmd
