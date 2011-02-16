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

#include <halmd/io/logger.hpp>
#include <halmd/observables/host/ssf.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using boost::fusion::at_key;
using namespace std;

namespace halmd
{
namespace observables { namespace host
{

template <int dimension>
ssf<dimension>::ssf(
    shared_ptr<density_modes_type> density_modes
  , unsigned int npart
)
  : _Base(density_modes)
  // dependency injection
  , density_modes(density_modes)
  // initialise members
  , npart_(npart)
{
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
 * compute SSF from sample of density Fourier modes
 */
template <int dimension>
void ssf<dimension>::compute_()
{
    scoped_timer<timer> timer_(at_key<sample_>(runtime_));

    typedef typename density_modes_type::density_modes_sample_type::mode_vector_type rho_vector_type;
    typedef typename density_modes_type::wavevectors_type::map_type wavevectors_map_type;
    typedef typename rho_vector_type::const_iterator rho_iterator;
    typedef typename wavevectors_map_type::const_iterator wavevector_iterator;
    typedef std::vector<accumulator<double> >::iterator result_iterator;

    // perform computation of partial SSF for all combinations of particle types
    wavevectors_map_type const& wavevectors = density_modes->wavevectors().values();
    unsigned int ntype = density_modes->ntype();
    unsigned int k = 0;
    for (unsigned char i = 0; i < ntype; ++i) {
        for (unsigned char j = i; j < ntype; ++j, ++k) {
            rho_iterator rho_q0 = density_modes->rho_sample->rho[i]->begin();
            rho_iterator rho_q1 = density_modes->rho_sample->rho[j]->begin();
            result_iterator result = result_accumulator_[k].begin();

            // iterate over ranges of wavevectors with equal magnitude
            wavevector_iterator q_begin = wavevectors.begin();
            while (q_begin != wavevectors.end()) {
                // find range with equal map key
                typedef typename wavevectors_map_type::value_type value_type;
                wavevector_iterator q_end = adjacent_find(
                    q_begin, wavevectors.end()
                  , bind(less<double>(), bind(&value_type::first, _1), bind(&value_type::first, _2))
                );
                // adjacent_find(begin, end, pred) returns the first iterator 'it' where 'pred(it, ++it)' holds
                // and 'end' if there is no match. Thus, the range of wavevectors in case
                // of a match is (q_begin, ++q_end). This is inconvenient and we fix it:
                if (q_end != wavevectors.end()) {
                    ++q_end;
                }

                // accumulate products of density modes with equal wavenumber
                double sum = 0;
                unsigned int count = 0;
                for (wavevector_iterator q = q_begin; q != q_end; ++q, ++rho_q0, ++rho_q1, ++count) {
                    // rho_q × rho_q^*
                    sum += real(*rho_q0) * real(*rho_q1) + imag(*rho_q0) * imag(*rho_q1);
                }
                // add result to accumulator
                (*result++)(sum / count);
                // start next range at end of current one
                q_begin = q_end;
            }
        }
    }
}

template <int dimension>
void ssf<dimension>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static string class_name("ssf_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                namespace_("host")
                [
                    class_<ssf, shared_ptr<_Base_Base>, bases<_Base_Base, _Base> >(class_name.c_str())
                        .def(constructor<
                            shared_ptr<density_modes_type>
                          , unsigned int
                        >())
                        .def("register_runtimes", &ssf::register_runtimes)
                ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
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

}} // namespace observables::host

} // namespace halmd
