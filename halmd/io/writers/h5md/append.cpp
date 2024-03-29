/*
 * Copyright © 2013-2014 Felix Höfling
 * Copyright © 2011      Peter Colberg
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
#include <boost/type_traits/has_dereference.hpp>
#include <luaponte/luaponte.hpp>
#include <luaponte/out_value_policy.hpp>
#include <memory>
#include <limits>
#include <stdexcept>
#include <stdint.h> // uint32_t, uint64_t
#include <type_traits>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/io/writers/h5md/append.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/raw_array.hpp>

using boost::multi_array;
using namespace std;

namespace halmd {
namespace io {
namespace writers {
namespace h5md {

append::append(
    H5::Group const& root
  , vector<string> const& location
  , std::shared_ptr<clock_type const> clock
)
  : clock_(clock)
  , last_step_(numeric_limits<int64_t>::lowest())
  , last_time_(numeric_limits<time_type>::lowest())
{
    if (location.size() < 1) {
        throw invalid_argument("group location");
    }
    group_ = h5xx::open_group(root, boost::join(location, "/"));
    step_dataset_ = h5xx::create_chunked_dataset<step_type>(group_, "step");
    time_dataset_ = h5xx::create_chunked_dataset<time_type>(group_, "time");
    group_.unlink("step");
    group_.unlink("time");
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , T const&
)
{
    return h5xx::create_chunked_dataset<T>(group, name);
}

template <typename T, typename Alloc>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , vector<T, Alloc> const& data
)
{
    return h5xx::create_chunked_dataset<vector<T, Alloc> >(group, name, data.size());
}

template <typename T>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , raw_array<T> const& data
)
{
    return h5xx::create_chunked_dataset<raw_array<T> >(group, name, data.size());
}

template <typename T, size_t N, typename Alloc>
static H5::DataSet create_dataset(
    H5::Group const& group
  , string const& name
  , multi_array<T, N, Alloc> const& data
)
{
    return h5xx::create_chunked_dataset<multi_array<T, N, Alloc> >(group, name, data.shape());
}

template <typename T>
static typename std::enable_if<boost::has_dereference<T>::type::value, void>::type
write_dataset(
    H5::DataSet& dataset
  , H5::Group const& group
  , string const& name
  , std::function<T ()> const& slot
)
{
    auto const& data = *slot();
    if (!h5xx::is_valid(dataset.getId())) {
        dataset = create_dataset(group, name, data);
    }
    h5xx::write_chunked_dataset(dataset, data);
}

template <typename T>
static typename std::enable_if<!boost::has_dereference<T>::type::value, void>::type
write_dataset(
    H5::DataSet& dataset
  , H5::Group const& group
  , string const& name
  , std::function<T ()> const& slot
)
{
    T data = slot();
    if (!h5xx::is_valid(dataset.getId())) {
        dataset = create_dataset(group, name, data);
    }
    h5xx::write_chunked_dataset(dataset, data);
}

template <typename T>
connection append::on_write(
    subgroup_type& group
  , std::function<T ()> const& slot
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    group = h5xx::open_group(group_, boost::join(location, "/"));
    h5xx::link(step_dataset_, group, "step");
    h5xx::link(time_dataset_, group, "time");
    return on_write_.connect(bind(&write_dataset<T>, H5::DataSet(), group, "value", slot));
}

template <typename T>
connection append::on_write_averaged(
    subgroup_type& group
  , std::function<T ()> const& value_slot
  , std::function<T ()> const& error_slot
  , std::function<uint64_t ()> const& count_slot
  , vector<string> const& location
)
{
    if (location.size() < 1) {
        throw invalid_argument("dataset location");
    }
    group = h5xx::open_group(group_, boost::join(location, "/"));
    h5xx::link(step_dataset_, group, "step");
    h5xx::link(time_dataset_, group, "time");

    H5::DataSet value_dataset, error_dataset, count_dataset;
    return on_write_.connect( [=]() mutable {
        write_dataset(value_dataset, group, "value", value_slot);
        write_dataset(error_dataset, group, "error", error_slot);
        write_dataset(count_dataset, group, "count", count_slot);
    });
}

connection append::on_prepend_write(slot_function_type const& slot)
{
    return on_prepend_write_.connect(slot);
}

connection append::on_append_write(slot_function_type const& slot)
{
    return on_append_write_.connect(slot);
}

void append::write()
{
    on_prepend_write_();
    write_step_time();
    on_write_();
    on_append_write_();
}

void append::write_step_time()
{
    step_type step = clock_->step();
    time_type time = clock_->time();
    if (static_cast<int64_t>(step) <= last_step_ || time < last_time_) {
        throw std::logic_error("Writing to H5MD file failed at step " + boost::lexical_cast<std::string>(step) +
                               "\nH5MD enforces a strictly increasing order.");
    }

    h5xx::write_chunked_dataset(step_dataset_, step);
    h5xx::write_chunked_dataset(time_dataset_, time);
    last_step_ = step;
    last_time_ = time;
}

static append::slot_function_type
wrap_write(std::shared_ptr<append> self)
{
    return [=]() {
        self->write();
    };
}

/**
 * As write slots we support functors with return by copy, return by const
 * reference and return by non-const reference. Non-const references are
 * supported since it is not possible to bind more than one data slot
 * under the same property name to Lua, therefore for modules with
 * read-write data (e.g. samples::host::phase_space) we bind only the
 * read-write slot returning a non-const reference.
 */
void append::luaopen(lua_State* L)
{
    using namespace boost::placeholders;
    using namespace luaponte;

    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("writers")
            [
                namespace_("h5md")
                [
                    class_<append, std::shared_ptr<append> >("append")
                        .def(constructor<H5::Group const&, vector<string> const&, std::shared_ptr<clock_type const> >())
                        .property("group", &append::group)
                        .property("write", &wrap_write)
                        .def("on_write", &append::on_write<float>, pure_out_value(_2))
                        .def("on_write", &append::on_write<float&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<float const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<double>, pure_out_value(_2))
                        .def("on_write", &append::on_write<double&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<double const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<uint32_t>, pure_out_value(_2))
                        .def("on_write", &append::on_write<uint32_t&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<uint32_t const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<uint64_t>, pure_out_value(_2))
                        .def("on_write", &append::on_write<uint64_t&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<uint64_t const&>, pure_out_value(_2))

                        .def("on_write", &append::on_write<fixed_vector<float, 2> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 2>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 2> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 3> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 3>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 3> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 6> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 6>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<float, 6> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 2> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 2>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 2> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 3> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 3>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 3> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 6> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 6>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<fixed_vector<double, 6> const&>, pure_out_value(_2))

                        .def("on_write", &append::on_write<vector<float> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<float>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<float> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<double> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<double>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<double> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<unsigned int> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<unsigned int>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<unsigned int> const&>, pure_out_value(_2))

                        .def("on_write", &append::on_write<vector<fixed_vector<float, 2> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 2> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 2> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 3> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 3> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 6> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 6> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<float, 6> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 2> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 2> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 2> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 3> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 3> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 6> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 6> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<fixed_vector<double, 6> > const&>, pure_out_value(_2))

                        // phase_space: position, velocity
                        .def("on_write", &append::on_write<raw_array<fixed_vector<float, 2>>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<float, 2>> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<float, 3>>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<float, 3>> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<float, 6>>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<float, 6>> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<double, 2>>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<double, 2>> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<double, 3>>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<double, 3>> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<double, 6>>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<fixed_vector<double, 6>> const&>, pure_out_value(_2))
                        // density_mode
                        .def("on_write", &append::on_write<std::shared_ptr<raw_array<fixed_vector<double, 2>> const>>, pure_out_value(_2))
                        .def("on_write", &append::on_write<std::shared_ptr<raw_array<fixed_vector<double, 3>> const>>, pure_out_value(_2))

                        .def("on_write", &append::on_write<raw_array<float>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<float> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<double>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<double> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<unsigned int>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<raw_array<unsigned int> const&>, pure_out_value(_2))

                        // ssf FIXME support output of accumulators as dataset triple (value, error, count)
                        .def("on_write", &append::on_write<raw_array<boost::array<double, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<float, 3> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<float, 3> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<float, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<float, 6> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<float, 6> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<float, 6> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<double, 3> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<double, 3> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<double, 3> > const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<double, 6> > >, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<double, 6> >&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<vector<boost::array<double, 6> > const&>, pure_out_value(_2))

                        .def("on_write", &append::on_write<multi_array<float, 2> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 2>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 2> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 3> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 3>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 3> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 6> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 6>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<float, 6> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 2> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 2>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 2> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 3> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 3>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 3> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 6> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 6>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<double, 6> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 2> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 2>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 2> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 3> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 3>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 3> const&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 6> >, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 6>&>, pure_out_value(_2))
                        .def("on_write", &append::on_write<multi_array<uint64_t, 6> const&>, pure_out_value(_2))

                        .def("on_write", &append::on_write_averaged<float>, pure_out_value(_2))
                        .def("on_write", &append::on_write_averaged<double>, pure_out_value(_2))

                        .def("on_write", &append::on_write_averaged<fixed_vector<float, 2>>, pure_out_value(_2))
                        .def("on_write", &append::on_write_averaged<fixed_vector<float, 3>>, pure_out_value(_2))
                        .def("on_write", &append::on_write_averaged<fixed_vector<float, 6>>, pure_out_value(_2))
                        .def("on_write", &append::on_write_averaged<fixed_vector<double, 2>>, pure_out_value(_2))
                        .def("on_write", &append::on_write_averaged<fixed_vector<double, 3>>, pure_out_value(_2))
                        .def("on_write", &append::on_write_averaged<fixed_vector<double, 6>>, pure_out_value(_2))

                        .def("on_prepend_write", &append::on_prepend_write)
                        .def("on_append_write", &append::on_append_write)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_writers_h5md_append(lua_State* L)
{
    append::luaopen(L);
    return 0;
}

} // namespace h5md
} // namespace writers
} // namespace io
} // namespace halmd
