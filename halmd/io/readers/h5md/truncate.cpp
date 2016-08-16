/*
 * Copyright © 2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <boost/algorithm/string/join.hpp> // boost::join
#include <luaponte/luaponte.hpp>
#include <luaponte/out_value_policy.hpp>
#include <stdexcept>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/readers/h5md/truncate.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace io {
namespace readers {
namespace h5md {

truncate::truncate(
    H5::Group const& root
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("group location");
    }
    group_ = root.openGroup(boost::join(location, "/"));
}

template <typename T>
connection truncate::on_read(
    subgroup_type& dataset
  , std::function<T ()> const& slot
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    dataset = group_.openDataSet(boost::join(location, "/"));
    return on_read_.connect(bind(&read_dataset<T>, dataset, slot));
}

connection truncate::on_prepend_read(slot_function_type const& slot)
{
    return on_prepend_read_.connect(slot);
}

connection truncate::on_append_read(slot_function_type const& slot)
{
    return on_append_read_.connect(slot);
}

void truncate::read()
{
    on_prepend_read_();
    on_read_();
    on_append_read_();
}

template <typename T>
void truncate::read_dataset(
    H5::DataSet dataset
  , std::function<T ()> const& slot
)
{
    h5xx::read_dataset(dataset, slot());
}

template<typename T>
void wrap_on_read(std::shared_ptr<truncate> self, H5::DataSet &group, std::function<void(const T&)> slot
        , std::vector<std::string> const& location)
{
    std::shared_ptr<T> array = std::make_shared<T>();

    self->on_read<T&>(group, [array]() -> T& { return *array; }, location);
    self->on_append_read([array, slot] () mutable {
        slot(*array);
        array.reset();
    });
}

void truncate::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("readers")
            [
                namespace_("h5md")
                [
                    class_<truncate, std::shared_ptr<truncate> >("truncate")
                        .def(constructor<H5::Group const&, vector<string> const&>())
                        .property("group", &truncate::group)
                        .def("read", &truncate::read)
                        .def("on_read", &truncate::on_read<float&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<double&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<fixed_vector<float, 2>&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<fixed_vector<float, 3>&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<fixed_vector<double, 2>&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<fixed_vector<double, 3>&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<float>&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<double>&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<unsigned int>&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<fixed_vector<float, 2> >&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<fixed_vector<float, 3> >&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<fixed_vector<double, 2> >&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<fixed_vector<double, 3> >&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<boost::array<float, 3> >&>, pure_out_value(_2))
                        .def("on_read", &truncate::on_read<vector<boost::array<double, 3> >&>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<float>>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<double>>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<unsigned int>>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<float, 2> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<float, 3> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<double, 2> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<double, 3> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<boost::array<float, 3> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<boost::array<double, 3> >>, pure_out_value(_2))
                        .def("on_prepend_read", &truncate::on_prepend_read)
                        .def("on_append_read", &truncate::on_append_read)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_readers_h5md_truncate(lua_State* L)
{
    truncate::luaopen(L);
    return 0;
}

} // namespace h5md
} // namespace readers
} // namespace io
} // namespace halmd
