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

#include <boost/bind.hpp>

#include <halmd/observables/observable.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables
{

template <typename observable_type>
typename signal<void ()>::slot_function_type
prepare_wrapper(shared_ptr<observable_type> observable)
{
    return bind(&observable_type::prepare, observable);
}

template <typename observable_type>
typename signal<void (double)>::slot_function_type
sample_wrapper(shared_ptr<observable_type> observable)
{
    return bind(&observable_type::sample, observable, _1);
}

template <int dimension>
void observable<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("observable_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            class_<observable, shared_ptr<observable> >(class_name.c_str())
                .def("register_observables", &observable::register_observables)
                .property("prepare", &prepare_wrapper<observable>)
                .property("sample", &sample_wrapper<observable>)
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &observable<3>::luaopen
    ]
    [
        &observable<2>::luaopen
    ];
}

} // namespace

} // namespace observables

} // namespace halmd
