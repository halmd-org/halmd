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

#include <boost/bind.hpp>
#include <functional>
#include <iterator>
#include <string>

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
    shared_ptr<density_modes_type> density_modes
  , unsigned int npart
)
  // dependency injection
  : density_modes(density_modes)
  // initialise members
  , npart_(npart)
  , time_(-1)
{
    // allocate memory
    unsigned int nq = density_modes->wavenumbers().size();
    unsigned int ntype = density_modes->value().size();
    unsigned int nssf = ntype * (ntype + 1) / 2; //< number of partial structure factors

    value_.resize(nssf);
    result_accumulator_.resize(nssf);
    for (unsigned int i = 0; i < nssf; ++i) {
        value_[i].resize(nq);
        result_accumulator_[i].resize(nq);
    }
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
    string root("structure/ssf/");
    // write wavenumbers only once
    writer.write_dataset(root + "wavenumbers", density_modes->wavenumbers(), "wavenumber grid");

    // register output writers for all partial structure factors
    unsigned char ntype = static_cast<unsigned char>(density_modes->value().size());
    assert('A' + density_modes->value().size() <= 'Z' + 1);
    unsigned int k = 0;
    for (unsigned char i = 0; i < ntype; ++i) {
        for (unsigned char j = i; j < ntype; ++j, ++k) {
            string label;
            label += 'A' + i;
            label += 'A' + j;
            writer.register_observable(
                root + label, &value_[k]
              , "partial static structure factor S_" + label + " (value, error, count)"
            );
        }
    }
}

/**
 * compute SSF from sample of density Fourier modes
 */
template <int dimension>
void ssf<dimension>::sample(double time)
{
    // acquire sample of density modes and compute SSF
    density_modes->acquire(time);
    compute_();

    // iterate over combinations of particle types
    for (unsigned int i = 0; i < value_.size(); ++i) {
        // transform accumulators to tuples (mean, error_of_mean, count)
        for (unsigned int j = 0; j < result_accumulator_[i].size(); ++j) {
            accumulator<double> const& acc = result_accumulator_[i][j];
            result_type& v = value_[i][j];
            v[0] = mean(acc);
            v[1] = (count(acc) > 1) ? error_of_mean(acc) : 0;
            v[2] = static_cast<double>(count(acc));
        }
    }
    time_ = time;   // store time for writer functions
}

/**
 * compute SSF from sample of density Fourier modes
 */
template <int dimension>
void ssf<dimension>::compute_()
{
    scoped_timer<timer> timer_(at_key<sample_>(runtime_));

    typedef typename density_modes_type::result_type::value_type::element_type rho_vector_type;
    typedef typename density_modes_type::wavevectors_type::map_type wavevectors_map_type;
    typedef typename rho_vector_type::const_iterator rho_iterator;
    typedef typename wavevectors_map_type::const_iterator wavevector_iterator;
    typedef std::vector<accumulator<double> >::iterator result_iterator;

    // perform computation of partial SSF for all combinations of particle types
    wavevectors_map_type const& wavevectors = density_modes->wavevectors().values();
    if (wavevectors.empty()) return; // nothing to do

    unsigned int ntype = density_modes->value().size();
    unsigned int k = 0;
    for (unsigned char i = 0; i < ntype; ++i) {
        for (unsigned char j = i; j < ntype; ++j, ++k) {
            rho_iterator rho_q0 = density_modes->value()[i]->begin();
            rho_iterator rho_q1 = density_modes->value()[j]->begin();
            result_iterator result = result_accumulator_[k].begin();

            // accumulate products of density modes with equal wavenumber,
            // iterate over sorted list of wavevectors
            double sum = 0;
            unsigned int count = 0;
            wavevector_iterator q = wavevectors.begin();
            wavevector_iterator q_next = q; ++q_next;
            while (q != wavevectors.end()) {
                // rho_q × rho_q^*
                sum += real(*rho_q0) * real(*rho_q1) + imag(*rho_q0) * imag(*rho_q1);
                ++count;
                // find end of range with equal wavenumber
                if (q_next == wavevectors.end() || q->first != q_next->first) {
                    (*result++)(sum / count);   // add result to output accumulator
                    sum = 0; count = 0;
                }
                ++q; ++q_next; ++rho_q0; ++rho_q1;
            }
        }
    }
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
                    .def(constructor<
                        shared_ptr<density_modes_type>
                      , unsigned int
                    >())
                    .def("register_runtimes", &ssf::register_runtimes)
                    .property("value", &ssf::value)
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &ssf<3>::luaopen
    ]
    [
        &ssf<2>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class ssf<3>;
template class ssf<2>;

} // namespace observables

} // namespace halmd
