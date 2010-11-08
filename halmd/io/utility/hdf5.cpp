/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <h5xx/h5xx.hpp>

#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{

static void luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "h5xx")
    [
        class_<H5::IdComponent>("id")

      , class_<H5::CommonFG>("common_fg")

      , class_<H5::AbstractDs>("abstract_dataset")

      , class_<H5::H5Object, H5::IdComponent>("object")

      , class_<H5::H5File, bases<H5::IdComponent, H5::CommonFG> >("file")
            .def("open_group", (H5::Group (*)(H5::CommonFG const&, string const&)) &h5xx::open_group)

      , class_<H5::Group, bases<H5::H5Object, H5::CommonFG> >("group")
            .def("open_group", (H5::Group (*)(H5::CommonFG const&, string const&)) &h5xx::open_group)

      , class_<H5::DataSet, bases<H5::H5Object, H5::AbstractDs> >("dataset")
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &luaopen
    ];
}

} // namespace halmd
