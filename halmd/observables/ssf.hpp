/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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
#include <lua.hpp>
#include <memory>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {

/**
 * Compute static structure factor for isotropic system.
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
    // typedef raw_array<accumulator<double>> result_type; FIXME
    typedef raw_array<boost::array<double, 3>> result_type;
    typedef std::shared_ptr<raw_array<fixed_vector<double, 2>> const> mode_type;
    typedef std::function<mode_type ()> mode_acquisitor_type;
    typedef observables::utility::wavevector<dimension> wavevector_type;
    typedef logger logger_type;

    /**
     * Construct static structure factor instance.
     *
     * The argument 'norm' is the normalisation factor (e.g. total number of particles).
     */
    ssf(
        mode_acquisitor_type mode1
      , mode_acquisitor_type mode2
      , std::shared_ptr<wavevector_type const> wavevector
      , double norm
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>("ssf")
    );

    /**
     * Compute SSF from samples of density Fourier modes.
     */
    result_type const& sample();

    /**
     * Functor wrapping sample() for class instance stored within std::shared_ptr
     */
    static std::function<result_type const& ()> sampler(std::shared_ptr<ssf> self)
    {
        return [=]() -> result_type const& {
            return self->sample();
        };
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** acquisitors yielding the density modes */
    mode_acquisitor_type mode1_;
    mode_acquisitor_type mode2_;
    /** wavevector grid */
    std::shared_ptr<wavevector_type const> wavevector_;
    /** logger instance */
    std::shared_ptr<logger_type> logger_;

    /** normalisation factor */
    double norm_;
    /** cached result for the static structure factor */
    result_type result_;

    /** observers for input data cached within std::shared_ptr */
    std::weak_ptr<typename mode_type::element_type> mode1_observer_;
    std::weak_ptr<typename mode_type::element_type> mode2_observer_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

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
