/*
 * Copyright © 2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_DYNAMICS_BLOCKING_SCHEME_HPP
#define HALMD_OBSERVABLES_DYNAMICS_BLOCKING_SCHEME_HPP

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

/**
 * Store input samples (phase space, density modes, ...) in
 * coarse-grained block structures and trigger computation of
 * time correlation functions.
 *
 * The module is driven by connecting to the signal on_sample,
 * time correlation functions are registered by the method add().
 */
class blocking_scheme
{
private:
    typedef halmd::signal<void ()> signal_type;

public:
    typedef samples::blocking_scheme_base block_sample_type;
    typedef signal_type::slot_function_type slot_function_type;
    typedef mdsim::clock clock_type;
    typedef clock_type::step_type step_type;
    typedef clock_type::time_type time_type;
    typedef logger logger_type;
    typedef boost::multi_array<time_type, 2> block_time_type;

    /**
     *  @param maximum_lag_time   maximum lag time for dynamic correlations
     *  @param resolution         time resolution of lowest level
     *  @param block_size         size of each block, determines coarse-graining factor
     *  @param shift              coarse-graining shift between odd and even levels,
     *                            if 0 it is computed as sqrt(block_size)
     */
    blocking_scheme(
        boost::shared_ptr<clock_type const> clock
      , double maximum_lag_time
      , double resolution
      , unsigned int block_size
      , unsigned int shift = 0
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /** add a time correlation function */
    connection on_correlate(boost::shared_ptr<correlation_base> tcf);
    /** add blocked input data, e.g., phase space points or density modes */
    connection on_sample(boost::shared_ptr<block_sample_type> block_sample);

    /** acquire and store current data from input sample,
     *  emit on_correlate_block for each full coarse-graining level
     */
    void sample();

    /** compute remaining correlations of partially filled coarse-graining levels */
    void finalise();

    /** signal triggers preparation of input sample */
    connection on_prepend_sample(slot_function_type const& slot);
    /** signal emitted after appending sample(s) to block */
    connection on_append_sample(slot_function_type const& slot);
    /** signal emitted before finalise */
    connection on_prepend_finalise(slot_function_type const& slot);
    /** signal emitted after finalise */
    connection on_append_finalise(slot_function_type const& slot);

    /** returns block size, i.e., the length of the coarse-graining levels */
    unsigned int block_size() const
    {
        return block_size_;
    }

    /** returns block count, i.e., the number of coarse-graining levels */
    unsigned int count() const
    {
        return interval_.size();
    }

    /** returns block-wise time grid for correlation functions */
    block_time_type const& time() const
    {
        return time_;
    }

    /** Lua bindings */
    static void luaopen(lua_State* L);

private:
    /** compute correlations and discard first entry */
    void process(unsigned int level);

     /** simulation clock */
    boost::shared_ptr<clock_type const> clock_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
    /** set of time correlation functions */
    slots<boost::shared_ptr<correlation_base> > tcf_;
    /** set of block structures holding the input samples (e.g., phase space point, density modes)
     *
     * We use std::set since it is a Unique Container.
     */
    slots<boost::shared_ptr<block_sample_type> > block_sample_;
    /** size (length) of the coarse-graining levels */
    unsigned int block_size_;
    /** sampling intervals for each coarse-graining level */
    std::vector<step_type> interval_;
    /** time grid of the resulting correlation functions */
    block_time_type time_;
    /** signal triggers preparation of input sample */
    signal_type on_prepend_sample_;
    /** signal emitted after appending sample(s) to block */
    signal_type on_append_sample_;
    /** signal emitted before finalise */
    signal_type on_prepend_finalise_;
    /** signal emitted after finalise */
    signal_type on_append_finalise_;
};

} // namespace dynamics
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_BLOCKING_SCHEME_HPP */
