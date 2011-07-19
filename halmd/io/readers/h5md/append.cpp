/*
 * Copyright © 2011  Peter Colberg
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

#include <boost/algorithm/string/join.hpp> // boost::join
#include <luabind/luabind.hpp>
#include <luabind/out_value_policy.hpp>
#include <stdexcept>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/readers/h5md/append.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace io {
namespace readers {
namespace h5md {

append::append(
    H5::Group const& root
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("group location");
    }
    group_ = root.openGroup(join(location, "/"));
}

template <typename T>
void append::on_read(
    subgroup_type& group
  , function<T ()> const& slot
  , vector<string> const& location
  , hsize_t index
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    group = h5xx::open_group(group_, join(location, "/"));
    H5::DataSet dataset = group.openDataSet("samples");
    on_read_.connect(bind(&read_dataset<T>, dataset, slot, index));
}

void append::on_prepend_read(signal<void (uint64_t)>::slot_function_type const& slot)
{
    on_prepend_read_.connect(slot);
}

void append::on_append_read(signal<void (uint64_t)>::slot_function_type const& slot)
{
    on_append_read_.connect(slot);
}

void append::read(uint64_t step)
{
    on_prepend_read_(step);
    on_read_();
    on_append_read_(step);
}

template <typename T>
void append::read_dataset(
    H5::DataSet dataset
  , function<T ()> const& slot
  , hsize_t index
)
{
    T data = slot();
    h5xx::read_dataset(dataset, &data, index);
}

static signal<void (uint64_t)>::slot_function_type
wrap_read(shared_ptr<append> instance)
{
    return bind(&append::read, instance, _1);
}

void append::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("readers")
            [
                namespace_("h5md")
                [
                    class_<append, shared_ptr<append> >("append")
                        .def(constructor<H5::Group const&, vector<string> const&>())
                        .property("read", &wrap_read)
                        .def("on_read", &append::on_read<float&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<double&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<float, 2>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<float, 3>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<double, 2>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<double, 3>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<float>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<double>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<float, 2> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<float, 3> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<double, 2> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<double, 3> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<array<float, 3> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<array<double, 3> >&>, pure_out_value(_2))
                        .def("on_prepend_read", &append::on_prepend_read)
                        .def("on_append_read", &append::on_append_read)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_readers_h5md_append(lua_State* L)
{
    append::luaopen(L);
    return 0;
}

} // namespace h5md
} // namespace readers
} // namespace io
} // namespace halmd
