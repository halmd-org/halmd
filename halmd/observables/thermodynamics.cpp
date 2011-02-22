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

#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace boost::fusion;
using namespace std;

namespace halmd
{
namespace observables
{

template <int dimension>
thermodynamics<dimension>::thermodynamics(
    shared_ptr<box_type> box
)
  // dependency injection
  : box(box)
{
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void thermodynamics<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * register observables
 */
template <int dimension>
void thermodynamics<dimension>::register_observables(writer_type& writer)
{
    writer.register_observable("TIME", &time_, "simulation time"); //< FIXME move time to mdsim::core
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
void thermodynamics<dimension>::sample(double time)
{
    scoped_timer<timer> timer_(at_key<sample_>(runtime_));
    en_pot_ = en_pot();
    en_kin_ = en_kin();
    v_cm_ = v_cm();
    en_tot_ = en_pot_ + en_kin_;
    temp_ = 2 * en_kin_ / dimension;
    density_ = box->density(); //< FIXME why is this duplicated in thermodynamics?
    pressure_ = density_ * (temp_ + virial() / dimension);
    hypervirial_ = hypervirial();
    time_ = time;
}

template <int dimension>
void thermodynamics<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("thermodynamics_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                class_<thermodynamics, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def("sample", &thermodynamics::sample)
                    .def("register_runtimes", &thermodynamics::register_runtimes)
                    .property("en_kin", &thermodynamics::en_kin)
                    .property("en_pot", &thermodynamics::en_pot)
                    .property("en_tot", &thermodynamics::en_tot)
                    .property("pressure", &thermodynamics::pressure)
                    .property("temp", &thermodynamics::temp)
                    .property("v_cm", &thermodynamics::v_cm)
                    .property("virial", &thermodynamics::virial)
                    .property("hypervirial", &thermodynamics::hypervirial)
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &thermodynamics<3>::luaopen
    ]
    [
        &thermodynamics<2>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class thermodynamics<3>;
template class thermodynamics<2>;

} // namespace observables

} // namespace halmd
