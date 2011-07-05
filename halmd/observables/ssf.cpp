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
#include <limits>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/observables/ssf.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {

template <int dimension>
ssf<dimension>::ssf(
    shared_ptr<density_mode_type> density_mode
  , unsigned int npart
)
  // dependency injection
  : density_mode(density_mode)
  // initialise members
  , npart_(npart)
  , step_(numeric_limits<uint64_t>::max())
{
    // allocate memory
    unsigned int nq = density_mode->wavenumber().size();
    unsigned int ntype = density_mode->value().size();
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
    profiler.register_runtime(runtime_.sample, "sample", "computation of static structure factor");
}

/**
 * register observables
 */
template <int dimension>
void ssf<dimension>::register_observables(writer_type& writer)
{
    string root("structure/ssf/");
    // write wavenumbers only once
    writer.write_dataset(root + "wavenumber", density_mode->wavenumber(), "wavenumber grid");

    // register output writers for all partial structure factors
    unsigned char ntype = static_cast<unsigned char>(density_mode->value().size());
    assert('A' + density_mode->value().size() <= 'Z' + 1);
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
void ssf<dimension>::sample(uint64_t step)
{
    if (step_ == step) {
        LOG_TRACE("[ssf] sample is up to date");
        return;
    }

    // acquire sample of density modes
    on_sample_(step);

    LOG_TRACE("[ssf] sampling");

    if (density_mode->step() != step) {
        throw logic_error("density modes sample was not updated");
    }

    // compute SSF
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
    step_ = step;   // store simulation step as time stamp
}

/**
 * compute SSF from sample of density Fourier modes
 */
template <int dimension>
void ssf<dimension>::compute_()
{
    scoped_timer<timer> timer_(runtime_.sample);

    typedef typename density_mode_type::result_type::value_type::element_type rho_vector_type;
    typedef typename density_mode_type::wavevector_type::map_type wavevector_map_type;
    typedef typename rho_vector_type::const_iterator rho_iterator;
    typedef typename wavevector_map_type::const_iterator wavevector_iterator;
    typedef std::vector<accumulator<double> >::iterator result_iterator;

    // perform computation of partial SSF for all combinations of particle types
    wavevector_map_type const& wavevector = density_mode->wavevector().value();
    if (wavevector.empty()) return; // nothing to do

    unsigned int ntype = density_mode->value().size();
    unsigned int k = 0;
    for (unsigned char i = 0; i < ntype; ++i) {
        for (unsigned char j = i; j < ntype; ++j, ++k) {
            rho_iterator rho_q0 = density_mode->value()[i]->begin();
            rho_iterator rho_q1 = density_mode->value()[j]->begin();
            result_iterator result = result_accumulator_[k].begin();

            // accumulate products of density modes with equal wavenumber,
            // iterate over sorted list of wavevectors
            wavevector_iterator q = wavevector.begin();
            wavevector_iterator q_next = q; ++q_next;
            while (q != wavevector.end()) {
                // compute Re[rho_q rho_q^*], add result to output accumulator
                double re = (real(*rho_q0) * real(*rho_q1) + imag(*rho_q0) * imag(*rho_q1));
                (*result)(re / npart_);
                // find end of range with equal wavenumber
                if (q_next == wavevector.end() || q->first != q_next->first) {
                    result++;   // next wavenumber: increment output iterator
                }
                ++q; ++q_next; ++rho_q0; ++rho_q1;
            }
        }
    }
}

template <typename ssf_type>
typename ssf_type::slot_function_type
sample_wrapper(shared_ptr<ssf_type> ssf)
{
    return bind(&ssf_type::sample, ssf, _1);
}

template <int dimension>
void ssf<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("ssf_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<ssf, shared_ptr<ssf> >(class_name.c_str())
                .def(constructor<
                    shared_ptr<density_mode_type>
                  , unsigned int
                >())
                .def("register_runtimes", &ssf::register_runtimes)
                .def("register_observables", &ssf::register_observables)
                .property("value", &ssf::value)
                .property("wavevector", &ssf::wavevector)
                .property("sample", &sample_wrapper<ssf>)
                .def("on_sample", &ssf::on_sample)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_ssf(lua_State* L)
{
    ssf<3>::luaopen(L);
    ssf<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class ssf<3>;
template class ssf<2>;

} // namespace observables
} // namespace halmd
