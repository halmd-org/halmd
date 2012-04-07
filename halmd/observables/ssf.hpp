/*
 * Copyright © 2011-2012  Felix Höfling and Peter Colberg
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

#ifndef HALMD_OBSERVABLES_SSF_HPP
#define HALMD_OBSERVABLES_SSF_HPP

#include <boost/array.hpp>
#include <boost/make_shared.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/samples/density_mode.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {

/**
 * Compute static structure factor.
 *
 * @f$ S_q^{(\alpha\beta)} = \langle \frac{1}{N} \rho_{\vec q}^{(\alpha)} \rho_{-\vec q}^{(\beta)} \rangle @f$
 * with the partial density modes
 * @f$ rho_{\vec q}^(\alpha) = \sum_{i=1}^{N_\alpha} \exp(i\vec q\cdot \vec r_i^{\alpha}) @f$
 * and the total number of particles @f$ N = \sum_\alpha N_\alpha @f$
 *
 * see e.g., Hansen & McDonald: Theory of simple liquids, chapter 4.1.
 */
template <int dimension>
class ssf
{
public:
    typedef observables::samples::density_mode<dimension> density_mode_type;
    typedef observables::utility::wavevector<dimension> wavevector_type;
    typedef mdsim::clock clock_type;
    typedef logger logger_type;
    typedef std::vector<boost::array<double, 3> > result_type;

    /**
     * Construct static structure factor instance.
     *
     * @param wavevector wavevector grid
     * @param norm normalisation factor (e.g. total number of particles)
     * @param clock simulation clock
     * @param logger logger instance
     */
    ssf(
        boost::shared_ptr<wavevector_type const> wavevector
      , double norm
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /**
     * Compute SSF from samples of density Fourier modes.
     */
    result_type const& sample(density_mode_type const& mode1, density_mode_type const& mode2);

    /**
     * Returns instance of wavevector class used to compute the SSF.
     */
    boost::shared_ptr<wavevector_type const> wavevector() const
    {
        return wavevector_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename clock_type::step_type step_type;
    typedef typename density_mode_type::mode_array_type rho_vector_type;

    /** wavevector grid */
    boost::shared_ptr<wavevector_type const> wavevector_;
    /** normalisation factor */
    double norm_;
    /** simulation clock */
    boost::shared_ptr<clock_type const> clock_;
    /** logger instance */
    boost::shared_ptr<logger_type> logger_;
    /** cached static structure factor */
    result_type result_;
    /** time stamp of data (simulation step) */
    step_type step_;

    typedef halmd::utility::profiler profiler_type;
    typedef profiler_type::accumulator_type accumulator_type;
    typedef profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type sample;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SSF_HPP */
