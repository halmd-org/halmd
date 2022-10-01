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

#include <halmd/config.hpp>

#include <algorithm>
#include <boost/algorithm/string/join.hpp> // boost::join
#include <cmath> // std::signbit
#include <limits>
#include <tuple>
#include <luaponte/luaponte.hpp>
#include <luaponte/out_value_policy.hpp>
#include <stdexcept>

#include <halmd/io/logger.hpp>
#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/readers/h5md/append.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/lua.hpp>

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
    group_ = root.openGroup(boost::join(location, "/"));
}

template <typename T>
connection append::on_read(
    subgroup_type& group
  , std::function<T ()> const& slot
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    group = h5xx::open_group(group_, boost::join(location, "/"));
    return on_read_.connect(bind(&read_dataset<T>, group, slot, boost::placeholders::_1));
}

connection append::on_prepend_read(slot_function_type const& slot)
{
    return on_prepend_read_.connect(slot);
}

connection append::on_append_read(slot_function_type const& slot)
{
    return on_append_read_.connect(slot);
}

void append::read_at_step(step_difference_type offset)
{
    on_prepend_read_();
    on_read_(bind(&read_step_index, offset, boost::placeholders::_1));
    on_append_read_();
}

void append::read_at_time(time_difference_type offset)
{
    on_prepend_read_();
    on_read_(bind(&read_time_index, offset, boost::placeholders::_1));
    on_append_read_();
}

template <typename T>
void append::read_dataset(
    H5::Group const& group
  , std::function<T ()> const& slot
  , index_function_type const& index
)
{
    H5::DataSet dataset = group.openDataSet("value");
    h5xx::read_chunked_dataset(dataset, slot(), index(group));
}

/**
 * Given a positive or negative step offset and a H5MD time series group,
 * this function returns the corresponding dataset index. If the offset is
 * negative, the last step incremented by one will be added to the offset.
 */
hsize_t append::read_step_index(
    step_difference_type offset
  , H5::Group const& group
)
{
    H5::DataSet dataset = group.openDataSet("step");
    std::vector<step_type> steps;
    h5xx::read_dataset(dataset, steps);
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

/**
 * Given a positive or negative time offset and a H5MD time series group,
 * this function returns the corresponding dataset index. If the offset is
 * negative, the last time will be added to the offset. As a special case,
 * note that the floating-point value 0.0 is different from -0.0, therefore
 * -0.0 may be used to refer exactly to the last time.
 *
 * When comparing two floating-point time values, we need to consider rounding
 * errors, as the total time calculated by mdsim::clock during the simulation
 * may be subject to other rounding errors than the time supplied by the caller
 * of append::read_at_time. To be strict but not too strict, we consider times
 * that differ less than a tolerance of 100 × (double precision floating-point
 * machine epsilon) × (minimum of two times) as equal.
 */
hsize_t append::read_time_index(
    time_difference_type offset
  , H5::Group const& group
)
{
    H5::DataSet dataset = group.openDataSet("time");
    std::vector<time_type> times;
    h5xx::read_dataset(dataset, times);
    if (times.size() < 1) {
        throw runtime_error("empty time dataset");
    }
    time_type time = signbit(offset) ? (offset + times.back()) : offset;
    std::vector<time_type>::const_iterator first, last;
    tie(first, last) = equal_range(
        times.begin()
      , times.end()
      , time
      , [](time_type time1, time_type time2) {
            return (time1 * (1 + 100 * numeric_limits<time_type>::epsilon()) < time2);
        }
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

/**
 * Wrapper function to allow one step reading.
 *
 * This wrapper internally creates an array and  connects it to the on_read slot.
 * Additionally it connects to the on_append_read slot and passes the read array
 * to the given slot function. After that the internal array is freed.
 * This way complete datasets can be directly passed to a single slot function
 * by const reference.
 */
template<typename T>
void wrap_on_read(
    std::shared_ptr<append> self
  , H5::Group& group
  , std::function<void(const T&)> slot
  , std::vector<std::string> const& location
)
{
    std::shared_ptr<T> array = std::make_shared<T>();

    self->on_read<T&>(group, [array]() -> T& { return *array; }, location);
    self->on_append_read([array, slot] () mutable {
        slot(*array);
        array.reset();
    });
}

void append::luaopen(lua_State* L)
{
    using namespace boost::placeholders;
    using namespace luaponte;

    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("readers")
            [
                namespace_("h5md")
                [
                    class_<append, std::shared_ptr<append> >("append")
                        .def(constructor<H5::Group const&, vector<string> const&>())
                        .property("group", &append::group)
                        .def("read_at_step", &append::read_at_step)
                        .def("read_at_time", &append::read_at_time)
                        .def("on_read", &append::on_read<float&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<double&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<float, 2>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<float, 3>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<double, 2>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<fixed_vector<double, 3>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<float>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<double>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<unsigned int>&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<float, 2> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<float, 3> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<double, 2> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<fixed_vector<double, 3> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<boost::array<float, 3> >&>, pure_out_value(_2))
                        .def("on_read", &append::on_read<vector<boost::array<double, 3> >&>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<float>>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<double>>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<unsigned int>>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<int>>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<float, 2> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<float, 3> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<double, 2> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<fixed_vector<double, 3> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<boost::array<float, 3> >>, pure_out_value(_2))
                        .def("on_read", &wrap_on_read<vector<boost::array<double, 3> >>, pure_out_value(_2))
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
