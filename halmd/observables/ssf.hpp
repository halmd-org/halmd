/*
 * Copyright © 2011-2012  Felix Höfling
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
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/density_mode.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {

/**
 * compute partial static structure factors
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
private:
    typedef halmd::signal<void ()> signal_type;

public:
    typedef observables::density_mode<dimension> density_mode_type;
    typedef mdsim::clock clock_type;
    typedef typename clock_type::step_type step_type;
    typedef logger logger_type;
    typedef typename signal_type::slot_function_type slot_function_type;

    typedef boost::array<double, 3> result_type;
    typedef fixed_vector<double, dimension> vector_type;

    static void luaopen(lua_State* L);

    ssf(
        std::vector<boost::shared_ptr<density_mode_type const> > const& density_mode
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    // compute ssf from sample of density Fourier modes and store with given simulation step
    void sample();

    connection on_sample(slot_function_type const& slot)
    {
        return on_sample_.connect(slot);
    }

    //! returns last computed values for static structure factor
    std::vector<std::vector<result_type> > const& value() const
    {
        return value_;
    }

    //! returns static structure factor for given type pair
    std::vector<result_type> const& value(unsigned int type1, unsigned int type2) const;

    //! returns instance of wavevector class used to compute the ssf
    typename density_mode_type::wavevector_type const& wavevector() const
    {
        return density_mode_.front()->wavevector();
    }

private:
    typedef halmd::utility::profiler profiler_type;
    typedef profiler_type::accumulator_type accumulator_type;
    typedef profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type sample;
    };

    std::vector<boost::shared_ptr<density_mode_type const> > density_mode_;
    boost::shared_ptr<clock_type const> clock_;
    boost::shared_ptr<logger_type> logger_;

    /** compute static structure factor and update result accumulators */
    void compute_();

    /** total number of particles, required for normalisation */
    unsigned int npart_;

    /**
     *  result for (partial) static structure factors
     *  in the order AA, AB, AC, …, BB, BC, …,CC
     *
     *  value_[i][j][0]:   mean value S_ab(k) for k = wavenumber_[j]
     *  value_[i][j][1]:   estimated error of mean
     *  value_[i][j][2]:   value count for the average
     */
    std::vector<std::vector<result_type> > value_;
    /** result accumulators */
    std::vector<std::vector<accumulator<double> > > result_accumulator_;
    /** time stamp of data (simulation step) */
    step_type step_;
    /** profiling runtime accumulators */
    runtime runtime_;

    signal_type on_sample_;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SSF_HPP */
