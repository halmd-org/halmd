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
#include <set>

#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd
{
namespace observables { namespace dynamics
{

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
public:
    typedef samples::blocking_scheme_base block_sample_type;
    typedef halmd::signal<void (uint64_t)> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    /**
     *  @param block_count  number of blocks, i.e., coarse-graining levels
     *  @param block_size   size of each block, determines coarse-graining factor
     *  @param shift        coarse-graining shift between odd and even levels
     *  @param resolution   time resolution of level 0
     */
    blocking_scheme(
        unsigned int block_count
      , unsigned int block_size
      , unsigned int shift
      , double resolution
    );

    /** add a time correlation function */
    void add_correlation(boost::shared_ptr<correlation_base> tcf)
    {
        tcf_.insert(tcf);
    }

    /** add blocked input data, e.g., phase space points or density modes */
    void add_data(boost::shared_ptr<block_sample_type> block_sample)
    {
        block_sample_.insert(block_sample);
    }

    /** acquire and store current data from input sample,
     *  emit on_correlate_block for each full coarse-graining level
     */
    void sample(uint64_t step);

    /** compute remaining correlations of partially filled coarse-graining levels */
    void finalise(uint64_t step);

    /** signal triggers preparation of input sample */
    void on_sample(slot_function_type const& slot)
    {
        on_sample_.connect(slot);
    }

private:
    /** set of time correlation functions */
    std::set<boost::shared_ptr<correlation_base> > tcf_;
    /** set of block structures holding the input samples (e.g., phase space point, density modes)
     *
     * We use std::set since it is a Unique Container.
     */
    std::set<boost::shared_ptr<block_sample_type> > block_sample_;

    /** sampling intervals for each coarse-graining level */
    std::vector<unsigned int> interval_;
    /** time grid of the resulting correlation functions */
    boost::multi_array<double, 2> time_;

    signal_type on_sample_;
};

}} // namespace observables::dynamics

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_BLOCKING_SCHEME_HPP */
