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
#include <halmd/utility/lua.hpp>
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
 * Resolve module dependencies
 */
template <int dimension>
void core<dimension>::depends()
{
    modules::depends<_Self, force_type>::required();
    modules::depends<_Self, neighbour_type>::required();
    modules::depends<_Self, sort_type>::optional();
    modules::depends<_Self, integrator_type>::required();
    modules::depends<_Self, particle_type>::required();
    modules::depends<_Self, position_type>::required();
    modules::depends<_Self, velocity_type>::required();
    modules::depends<_Self, profiler_type>::required();
}

template <int dimension>
void core<dimension>::select(po::options const& vm)
{
    if (vm["dimension"].as<int>() != dimension) {
        throw unsuitable_module("mismatching option dimension");
    }
}

/**
 * Initialize simulation
 */
template <int dimension>
core<dimension>::core(modules::factory& factory, po::options const& vm)
  // dependency injection
  : profiler(modules::fetch<profiler_type>(factory, vm))
  , force(modules::fetch<force_type>(factory, vm))
  , neighbour(modules::fetch<neighbour_type>(factory, vm))
  , sort(modules::fetch<sort_type>(factory, vm))
  , integrator(modules::fetch<integrator_type>(factory, vm))
  , particle(modules::fetch<particle_type>(factory, vm))
  , position(modules::fetch<position_type>(factory, vm))
  , velocity(modules::fetch<velocity_type>(factory, vm))
  // initialise attributes
  , step_counter_(0)
{
    /*@{ FIXME remove pre-Lua hack */
    register_runtimes(*profiler);
    /*@}*/
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
static luabind::scope register_lua(char const* class_name)
{
    using namespace luabind;
    return
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<T, shared_ptr<T> >(class_name)
                    .scope
                    [
                        def("options", &T::options)
                    ]
            ]
        ];
}

static lua_registry::iterator dummy = (
    lua_registry::get()->push_back( register_lua<core<3> >("core_3_") )
  , lua_registry::get()->push_back( register_lua<core<2> >("core_2_") )
  , lua_registry::get()->end()
);

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

template class module<mdsim::core<3> >;
template class module<mdsim::core<2> >;

} // namespace halmd
