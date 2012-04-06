/*
 * Copyright © 2011-2012  Peter Colberg
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
#include <stdint.h> // uint64_t

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/writers/h5md/truncate.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/raw_allocator.hpp>

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
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , T const&
)
{
    return h5xx::create_dataset<T>(group, name);
}

template <typename T, typename Alloc>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , vector<T, Alloc> const& data
)
{
    return h5xx::create_dataset<vector<T, Alloc> >(group, name, data.size());
}

template <typename T, size_t N, typename Alloc>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , multi_array<T, N, Alloc> const& data
)
{
    return h5xx::create_dataset<multi_array<T, N, Alloc> >(group, name, data.shape());
}

template <typename T>
static void write_dataset(
    H5::DataSet& dataset
  , H5::Group const& group
  , string const& name
  , function<T ()> const& slot
)
{
    T const& data = slot();
    if (!dataset.getId()) {
        dataset = create_dataset(group, name, data);
    }
    h5xx::write_dataset(dataset, data);
}

template <typename T>
connection truncate::on_write(
    subgroup_type& dataset
  , function<T ()> const& slot
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    return on_write_.connect(bind(&write_dataset<T>, dataset, group_, join(location, "/"), slot));
}

connection truncate::on_prepend_write(slot_function_type const& slot)
{
    return on_prepend_write_.connect(slot);
}

connection truncate::on_append_write(slot_function_type const& slot)
{
    return on_append_write_.connect(slot);
}

void truncate::write()
{
    on_prepend_write_();
    on_write_();
    on_append_write_();
}

static truncate::slot_function_type
wrap_write(shared_ptr<truncate> instance)
{
    return bind(&truncate::write, instance);
}

/**
 * As write slots we support functors with return by copy, return by const
 * reference and return by non-const reference. Non-const references are
 * supported since it is not possible to bind more than one data slot
 * under the same property name to Lua, therefore for modules with
 * read-write data (e.g. samples::host::phase_space) we bind only the
 * read-write slot returning a non-const reference.
 */
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
                        .property("group", &truncate::group)
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
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 2>, raw_allocator<fixed_vector<float, 2> > >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 2>, raw_allocator<fixed_vector<float, 2> > > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 3>, raw_allocator<fixed_vector<float, 3> > >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<float, 3>, raw_allocator<fixed_vector<float, 3> > > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 2>, raw_allocator<fixed_vector<double, 2> > >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 2>, raw_allocator<fixed_vector<double, 2> > > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 3>, raw_allocator<fixed_vector<double, 3> > >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<fixed_vector<double, 3>, raw_allocator<fixed_vector<double, 3> > > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<unsigned int, raw_allocator<unsigned int > >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<unsigned int, raw_allocator<unsigned int > > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<float, 3> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<float, 3> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<float, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<double, 3> > >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<double, 3> >&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<vector<array<double, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<float, 2> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<float, 2>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<float, 2> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<float, 3> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<float, 3>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<float, 3> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<double, 2> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<double, 2>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<double, 2> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<double, 3> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<double, 3>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<double, 3> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<uint64_t, 2> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<uint64_t, 2>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<uint64_t, 2> const&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<uint64_t, 3> >, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<uint64_t, 3>&>, pure_out_value(_2))
                        .def("on_write", &truncate::on_write<multi_array<uint64_t, 3> const&>, pure_out_value(_2))
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
