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

#include <boost/algorithm/string/join.hpp>
#include <boost/range/iterator_range.hpp>
#include <luabind/luabind.hpp>
#include <luabind/out_value_policy.hpp>
#include <stdexcept>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/writers/h5md/truncate.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::algorithm; // join
using namespace std;

namespace halmd {
namespace io {
namespace writers {
namespace h5md {

truncate::truncate(
    file_type const& file
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("group location");
    }
    group_ = h5xx::open_group(file.root(), join(location, "/"));
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , function<T ()> const&
)
{
    return h5xx::create_dataset<T>(group, name, 1);
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , function<T const& ()> const&
)
{
    return h5xx::create_dataset<T>(group, name, 1);
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , function<T& ()> const&
)
{
    return h5xx::create_dataset<T>(group, name, 1);
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , function<vector<T> ()> const& slot
)
{
    vector<T> data = slot();
    return h5xx::create_dataset<vector<T> >(group, name, data.size(), 1);
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , function<vector<T> const& ()> const& slot
)
{
    vector<T> const& data = slot();
    return h5xx::create_dataset<vector<T> >(group, name, data.size(), 1);
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , function<vector<T>& ()> const& slot
)
{
    vector<T>& data = slot();
    return h5xx::create_dataset<vector<T> >(group, name, data.size(), 1);
}

template <typename slot_type>
static void write_dataset(H5::DataSet dataset, slot_type const& slot)
{
    typename slot_type::result_type data = slot();
    h5xx::write_dataset(dataset, data, 0);
}

template <typename slot_type>
void truncate::on_write(
    H5::DataSet& dataset
  , slot_type const& slot
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    H5::Group group;
    if (location.size() > 1) {
        group = h5xx::open_group(
            group_
          , join(make_iterator_range(location.begin(), location.end() - 1), "/")
        );
    }
    else {
        group = group_;
    }
    dataset = create_dataset(group, location.back(), slot);
    on_write_.connect(bind(&write_dataset<slot_type>, dataset, slot));
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
                        .def(constructor<file const&, vector<string> const&>())
                        .property("write", &wrap_write)
                        .def("on_write", &truncate::on_write<function<float ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<double ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<float const&()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<double const&()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<float, 2> ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<float, 3> ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<double, 2> ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<double, 3> ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<float, 2> const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<float, 3> const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<double, 2> const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<fixed_vector<double, 3> const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<float> ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<double> ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<float> const&()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<double> const&()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<float, 2> > ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<float, 3> > ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<double, 2> > ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<double, 3> > ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<float, 2> > const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<float, 3> > const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<double, 2> > const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<fixed_vector<double, 3> > const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<array<float, 3> > ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<array<double, 3> > ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<array<float, 3> > const& ()> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<function<vector<array<double, 3> > const& ()> >, pure_out_value(_2))
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
