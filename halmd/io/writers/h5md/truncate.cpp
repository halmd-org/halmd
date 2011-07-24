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
#include <halmd/io/writers/h5md/truncate.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace io {
namespace writers {
namespace h5md {

truncate::truncate(
    H5::Group const& root
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("group location");
    }
    group_ = h5xx::open_group(root, join(location, "/"));
}

template <typename T>
void truncate::on_write(
    subgroup_type& dataset
  , function<T ()> const& slot
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    dataset = create_dataset(group_, join(location, "/"), slot);
    on_write_.connect(bind(&write_dataset<T>, dataset, slot));
}

void truncate::on_prepend_write(signal<void (uint64_t)>::slot_function_type const& slot)
{
    on_prepend_write_.connect(slot);
}

void truncate::on_append_write(signal<void (uint64_t)>::slot_function_type const& slot)
{
    on_append_write_.connect(slot);
}

void truncate::write(uint64_t step)
{
    on_prepend_write_(step);
    on_write_();
    on_append_write_(step);
}

template <typename T>
H5::DataSet truncate::create_dataset(
    H5::Group const& group
  , string const& name
  , function<T ()> const&
)
{
    return h5xx::create_unique_dataset<T>(group, name);
}

template <typename T>
H5::DataSet truncate::create_dataset(
    H5::Group const& group
  , string const& name
  , function<T const& ()> const&
)
{
    return h5xx::create_unique_dataset<T>(group, name);
}

template <typename T>
H5::DataSet truncate::create_dataset(
    H5::Group const& group
  , string const& name
  , function<T& ()> const&
)
{
    return h5xx::create_unique_dataset<T>(group, name);
}

template <typename T>
H5::DataSet truncate::create_dataset(
    H5::Group const& group
  , string const& name
  , function<vector<T> ()> const& slot
)
{
    vector<T> data = slot();
    return h5xx::create_unique_dataset<vector<T> >(group, name, data.size());
}

template <typename T>
H5::DataSet truncate::create_dataset(
    H5::Group const& group
  , string const& name
  , function<vector<T> const& ()> const& slot
)
{
    vector<T> const& data = slot();
    return h5xx::create_unique_dataset<vector<T> >(group, name, data.size());
}

template <typename T>
H5::DataSet truncate::create_dataset(
    H5::Group const& group
  , string const& name
  , function<vector<T>& ()> const& slot
)
{
    vector<T>& data = slot();
    return h5xx::create_unique_dataset<vector<T> >(group, name, data.size());
}

template <typename T>
void truncate::write_dataset(
    H5::DataSet dataset
  , function<T ()> const& slot
)
{
    T data = slot();
    h5xx::write_unique_dataset(dataset, data);
}

static signal<void (uint64_t)>::slot_function_type
wrap_write(shared_ptr<truncate> instance)
{
    return bind(&truncate::write, instance, _1);
}

void truncate::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("writers")
            [
                namespace_("h5md")
                [
                    class_<truncate, shared_ptr<truncate> >("truncate")
                        .def(constructor<H5::Group const&, vector<string> const&>())
                        .property("write", &wrap_write)
                        .def("on_write", &truncate::on_write<float>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<float&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<float const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<double>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<double&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<double const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<float, 2> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<float, 2>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<float, 2> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<float, 3> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<float, 3>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<float, 3> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<double, 2> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<double, 2>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<double, 2> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<double, 3> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<double, 3>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<fixed_vector<double, 3> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<float> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<float>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<float> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<double> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<double>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<double> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 2> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 2> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 2> > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 3> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 3> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 2> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 2> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 2> > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 3> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 3> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<float, 3> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<float, 3> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<float, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<double, 3> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<double, 3> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<double, 3> > const&>, pure_out_value(_2))
                        .def("on_prepend_write", &truncate::on_prepend_write)
                        .def("on_append_write", &truncate::on_append_write)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_writers_h5md_truncate(lua_State* L)
{
    truncate::luaopen(L);
    return 0;
}

} // namespace h5md
} // namespace writers
} // namespace io
} // namespace halmd
