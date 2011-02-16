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

#ifndef HALMD_OBSERVABLES_HOST_SSF_HPP
#define HALMD_OBSERVABLES_HOST_SSF_HPP

#include <boost/tuple/tuple.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/host/density_modes.hpp>
#include <halmd/observables/ssf.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace observables { namespace host
{

/**
 * compute static structure factor
 */

template <int dimension>
class ssf
  : public observables::ssf<dimension>
{
public:
    typedef observables::ssf<dimension> _Base;
#ifndef USE_HOST_SINGLE_PRECISION
    typedef observables::host::density_modes<dimension, double> density_modes_type;
#else
    typedef observables::host::density_modes<dimension, float> density_modes_type;
#endif
    typedef halmd::utility::profiler profiler_type;
    typedef fixed_vector<double, dimension> vector_type;

    boost::shared_ptr<density_modes_type> density_modes;

    static void luaopen(lua_State* L);

    ssf(
        boost::shared_ptr<density_modes_type> density_modes
      , unsigned int npart
    );

    void register_runtimes(profiler_type& profiler);

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG( sample_, "computation of static structure factor" );

protected:
    /** compute static structure factor and update result accumulators */
    virtual void compute_();

    /** result accumulators */
    using _Base::result_accumulator_;

    /** total number of particles, required for normalisation */
    unsigned int npart_;

    // list of profiling timers
    boost::fusion::map<
        boost::fusion::pair<sample_, accumulator<double> >
    > runtime_;
};

}} // namespace observables::host

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_SSF_HPP */
