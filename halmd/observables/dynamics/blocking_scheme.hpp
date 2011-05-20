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

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd
{
namespace observables { namespace dynamics
{

/**
 * Store input samples (phase space, density modes, ...) in a
 * coarse-grained block structure and provide a signal
 * on_correlate_block which correlation functions connect to.
 */
class blocking_scheme
{
public:
    typedef samples::blocking_scheme_base block_sample_type;
    typedef halmd::signal<void (uint64_t)> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    /**
     *  @param block_sample shared pointer to block structure of input samples
     *  @param factor       coarse-graining factor between levels (i, i+2)
     *  @param shift        coarse-graining shift between odd and even levels
     *  @param resolution   time resolution of level 0
     */
    blocking_scheme(
        boost::shared_ptr<block_sample_type> block_sample
      , unsigned int factor
      , unsigned int shift
      , double resolution
    );

    /** acquire and store current data from input sample,
     *  emit on_correlate_block for each full coarse-graining level
     */
    void acquire(uint64_t step);

    /** compute remaining correlations of partially filled coarse-graining levels */
    void finalise(uint64_t);

    /** signal is emitted for preparation of input sample */
    void on_acquire(slot_function_type const& slot)
    {
        on_acquire_.connect(slot);
    }

    /** signal is emitted on a full block and carries its level index */
    void on_correlate_block(slot_function_type const& slot)
    {
        on_correlate_block_.connect(slot);
    }

private:
    boost::shared_ptr<block_sample_type> block_sample_;

    /** sampling intervals for each coarse-graining level */
    std::vector<unsigned int> interval_;
    /** time grid of the resulting correlation functions */
    boost::multi_array<double, 2> time_;

    signal_type on_acquire_;
    signal_type on_correlate_block_;
};

}} // namespace observables::dynamics

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_BLOCKING_SCHEME_HPP */
