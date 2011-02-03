/*
 * Copyright © 2011  Felix Höfling
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

#include <iterator>

#include <halmd/algorithm/host/pick_lattice_points.hpp>
#include <halmd/observables/ssf.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using boost::fusion::at_key;
using namespace std;

namespace halmd
{
namespace observables
{

template <int dimension>
ssf<dimension>::ssf(
    shared_ptr<box_type> box
  , std::vector<double> const& wavenumbers
  , double tolerance
  , unsigned int max_count
)
  // dependency injection
  : box(box)
  // initialise members
  , wavenumbers_(wavenumbers)
  , tolerance_(tolerance)
  , max_count_(max_count)
  , time_(-1)
{
    algorithm::host::pick_lattice_points_from_shell(
        wavenumbers.begin(), wavenumbers.end()
      , inserter(wavevectors_, wavevectors_.begin())
      , element_div(vector_type(2 * M_PI), box->length())
      , tolerance
      , max_count
    );

    // remove wavenumbers with no compatible wavevectors
    for (vector<double>::iterator q_it = wavenumbers_.begin(); q_it != wavenumbers_.end(); ++q_it) {
        if (!wavevectors_.count(*q_it)) {
            LOG_WARNING("No wavevector compatible with |q| ≈ " << *q_it << ". Value discarded");
            wavenumbers_.erase(q_it--);   // post-decrement iterator, increment at end of loop
        }
    }

    // allocate memory
    result_.resize(wavenumbers_.size());
    result_acc_.resize(wavenumbers_.size());
    // FIXME support HDF5 output of tuples
    value_.resize(wavenumbers_.size());
    error_.resize(wavenumbers_.size());
    count_.resize(wavenumbers_.size());
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void ssf<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * register observables
 */
template <int dimension>
void ssf<dimension>::register_observables(writer_type& writer)
{
    // write wavenumbers only once
    writer.write_dataset("structure/ssf/wavenumbers", wavenumbers_, "wavenumber grid");

    // register output writers
    // FIXME support HDF5 output of tuples
    // writer.register_observable_once("structure/ssf", &result_, "static structure factor (value, error, count)");
    writer.register_observable("structure/ssf/value", &value_, "static structure factor (mean value)");
    writer.register_observable("structure/ssf/error", &error_, "error of mean");
    writer.register_observable("structure/ssf/count", &count_, "count of averaged values");
}

/**
 * Sample all ssf
 */
template <int dimension>
void ssf<dimension>::sample(double time)
{
    scoped_timer<timer> timer_(at_key<sample_>(runtime_));
    compute_();

    // transform accumulators to tuples (mean, error_of_mean, count)
    for (unsigned int i = 0; i < result_acc_.size(); ++i) {
        accumulator<double> const& a = result_acc_[i];
        unsigned int c = count(a);
        result_[i] = make_tuple(mean(a), (c > 1) ? error_of_mean(a) : 0, c);
    }

    // FIXME support HDF5 output of tuples
    // split vector of tuples
    double const& (*get0)(tuple<double, double, unsigned>::inherited const&) = &boost::get<0>;
    double const& (*get1)(tuple<double, double, unsigned>::inherited const&) = &boost::get<1>;
    unsigned const& (*get2)(tuple<double, double, unsigned>::inherited const&) = &boost::get<2>;
    transform(result_.begin(), result_.end(), value_.begin(), bind(get0, _1));
    transform(result_.begin(), result_.end(), error_.begin(), bind(get1, _1));
    transform(result_.begin(), result_.end(), count_.begin(), bind(get2, _1));
    time_ = time;   // store time for writer functions
}

template <int dimension>
void ssf<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("ssf_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                class_<ssf, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def("sample", &ssf::sample)
                    .def("register_runtimes", &ssf::register_runtimes)
                    .property("wavenumbers", &ssf::wavenumbers)
                    .property("wavevectors", &ssf::wavevectors)
                    .property("result", &ssf::result)
                    .property("tolerance", &ssf::tolerance)
                    .property("maximum_count", &ssf::maximum_count)
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &ssf<3>::luaopen
    ]
    [
        &ssf<2>::luaopen
    ];
}

// explicit instantiation
template class ssf<3>;
template class ssf<2>;

} // namespace observables

} // namespace halmd
