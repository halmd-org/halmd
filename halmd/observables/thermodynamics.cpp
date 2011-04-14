/*
 * Copyright © 2010  Felix Höfling
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

#include <limits>

#include <halmd/io/logger.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables
{

template <int dimension>
thermodynamics<dimension>::thermodynamics(
    shared_ptr<box_type> box
  , shared_ptr<clock_type> clock
)
  // dependency injection
  : box(box)
  , clock(clock)
  // initialise members
  , step_(numeric_limits<uint64_t>::max())
{
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void thermodynamics<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.sample, "sample", "computation of macroscopic state variables");
}

/**
 * register observables
 */
template <int dimension>
void thermodynamics<dimension>::register_observables(writer_type& writer)
{
    writer.register_observable("TIME", &time_, "simulation time"); //< FIXME move time to mdsim::clock
    writer.register_observable("EPOT", &en_pot_, "mean potential energy per particle");
    writer.register_observable("EKIN", &en_kin_, "mean kinetic energy per particle");
    writer.register_observable("ETOT", &en_tot_, "mean total energy per particle");
    writer.register_observable("VCM", &v_cm_, "centre-of-mass velocity");
    writer.register_observable("PRESS", &pressure_, "virial pressure");
    writer.register_observable("XVIR", &hypervirial_, "hypervirial sum ");
    writer.register_observable("TEMP", &temp_, "temperature");
}

/**
 * Sample macroscopic state variables
 *
 * Compute state variables and take care that expensive functions are
 * called only once.
 */
template <int dimension>
void thermodynamics<dimension>::sample(uint64_t step)
{
    if (step_ == step) {
        LOG_TRACE("[thermodynamics] sample is up to date");
        return;
    }

    LOG_TRACE("[thermodynamics] acquire sample");

    scoped_timer<timer> timer_(runtime_.sample);
    en_pot_ = en_pot();
    en_kin_ = en_kin();
    v_cm_ = v_cm();
    en_tot_ = en_pot_ + en_kin_;
    temp_ = 2 * en_kin_ / dimension;
    density_ = box->density(); //< FIXME why is this duplicated in thermodynamics?
    pressure_ = density_ * (temp_ + virial() / dimension);
    hypervirial_ = hypervirial();
    time_ = clock->time();
    step_ = step;
}

template <typename thermodynamics_type>
typename thermodynamics_type::slot_function_type
prepare_wrapper(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::prepare, thermodynamics);
}

template <typename thermodynamics_type>
typename thermodynamics_type::slot_function_type
sample_wrapper(shared_ptr<thermodynamics_type> thermodynamics)
{
    return bind(&thermodynamics_type::sample, thermodynamics, _1);
}

template <int dimension>
void thermodynamics<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("thermodynamics_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<thermodynamics, shared_ptr<thermodynamics> >(class_name.c_str())
                .def("register_runtimes", &thermodynamics::register_runtimes)
                .def("register_observables", &thermodynamics::register_observables)
                .property("prepare", &prepare_wrapper<thermodynamics>)
                .property("sample", &sample_wrapper<thermodynamics>)
                .property("en_kin", &thermodynamics::en_kin)
                .property("en_pot", &thermodynamics::en_pot)
                .property("en_tot", &thermodynamics::en_tot)
                .property("pressure", &thermodynamics::pressure)
                .property("temp", &thermodynamics::temp)
                .property("v_cm", &thermodynamics::v_cm)
                .property("virial", &thermodynamics::virial)
                .property("hypervirial", &thermodynamics::hypervirial)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_thermodynamics(lua_State* L)
{
    thermodynamics<3>::luaopen(L);
    thermodynamics<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class thermodynamics<3>;
template class thermodynamics<2>;

} // namespace observables

} // namespace halmd
