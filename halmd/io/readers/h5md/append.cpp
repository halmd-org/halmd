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

#include <algorithm>
#include <boost/algorithm/string/join.hpp> // boost::join
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/tuple/tuple.hpp> // boost::tie
#include <cmath> // std::signbit
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
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    group = h5xx::open_group(group_, join(location, "/"));
    H5::DataSet dataset = group.openDataSet("samples");
    on_read_.connect(bind(&read_dataset<T>, dataset, slot, _1, group));
}

void append::on_prepend_read(slot_function_type const& slot)
{
    on_prepend_read_.connect(slot);
}

void append::on_append_read(slot_function_type const& slot)
{
    on_append_read_.connect(slot);
}

void append::read_step(step_difference_type offset)
{
    on_prepend_read_();
    on_read_(bind(&read_step_index, offset, _1));
    on_append_read_();
}

void append::read_time(time_difference_type offset)
{
    on_prepend_read_();
    on_read_(bind(&read_time_index, offset, _1));
    on_append_read_();
}

template <typename T>
void append::read_dataset(
    H5::DataSet dataset
  , function<T ()> const& slot
  , index_function_type const& index
  , H5::Group const& group
)
{
    T data = slot();
    h5xx::read_dataset(dataset, &data, index(group));
}

hsize_t append::read_step_index(
    step_difference_type offset
  , H5::Group const& group
)
{
    H5::DataSet dataset = group.openDataSet("step");
    std::vector<step_type> steps;
    h5xx::read_unique_dataset(dataset, &steps);
    if (steps.size() < 1) {
        throw runtime_error("empty step dataset");
    }
    step_type step = (offset < 0) ? (offset + steps.back() + 1) : offset;
    std::vector<step_type>::const_iterator first, last;
    tie(first, last) = equal_range(steps.begin(), steps.end(), step);
    if (first == last) {
        LOG_ERROR("no step " << step << " in dataset " << h5xx::path(dataset));
        throw domain_error("nonexistent step");
    }
    if (last - first > 1) {
        LOG_ERROR("ambiguous step " << step << " in dataset " << h5xx::path(dataset));
        throw domain_error("ambiguous step");
    }
    LOG("reading " << h5xx::path(group) << " at step " << *first);
    return first - steps.begin();
}

hsize_t append::read_time_index(
    time_difference_type offset
  , H5::Group const& group
)
{
    H5::DataSet dataset = group.openDataSet("time");
    time_type timestep = h5xx::read_attribute<time_type>(dataset, "timestep");
    std::vector<time_type> times;
    h5xx::read_unique_dataset(dataset, &times);
    if (times.size() < 1) {
        throw runtime_error("empty time dataset");
    }
    time_type time = signbit(offset) ? (offset + times.back()) : offset;
    std::vector<time_type>::const_iterator first, last;
    tie(first, last) = equal_range(
        times.begin()
      , times.end()
      , time
      , lambda::_1 + (timestep / 2) < lambda::_2
    );
    if (first == last) {
        LOG_ERROR("no time " << time << " in dataset " << h5xx::path(dataset));
        throw domain_error("nonexistent time");
    }
    if (last - first > 1) {
        LOG_ERROR("ambiguous time " << time << " in dataset " << h5xx::path(dataset));
        throw domain_error("ambiguous time");
    }
    LOG("reading " << h5xx::path(group) << " at time " << *first);
    return first - times.begin();
}

static append::slot_function_type
wrap_read_step(shared_ptr<append> instance, append::step_difference_type offset)
{
    return bind(&append::read_step, instance, offset);
}

static append::slot_function_type
wrap_read_time(shared_ptr<append> instance, append::time_difference_type offset)
{
    return bind(&append::read_time, instance, offset);
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
                        .def("read_step", &wrap_read_step)
                        .def("read_time", &wrap_read_time)
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
