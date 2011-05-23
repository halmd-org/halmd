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

#ifndef HALMD_OBSERVABLES_DYNAMICS_CORRELATION_HPP
#define HALMD_OBSERVABLES_DYNAMICS_CORRELATION_HPP

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

#include <halmd/observables/samples/blocking_scheme.hpp>

namespace halmd
{
namespace observables { namespace dynamics
{

/**
 * Store input samples (phase space, density modes, ...) in a
 * coarse-grained block structure and provide a signal
 * on_correlate_block which correlation functions connect to.
 */
class correlation_base
{
public:
    correlation_base() {}
    virtual ~correlation_base() {}

    /** compute correlations at the given coarse-graining level */
    virtual void compute(unsigned int level) = 0;
};

template <typename tcf_type>
class correlation
  : public correlation_base
{
    typedef typename tcf_type::sample_type sample_type;
    typedef typename tcf_type::result_type result_type;
    typedef observables::samples::blocking_scheme<sample_type> block_sample_type;

    correlation(boost::shared_ptr<block_sample_type> block_sample);
    virtual ~correlation() {}

    virtual void compute(unsigned int level);

private:
    boost::shared_ptr<block_sample_type> block_sample_;
    boost::multi_array<result_type, 2> result_;
};

template <typename tcf_type>
correlation<tcf_type>::correlation(boost::shared_ptr<block_sample_type> block_sample)
  // dependency injection
  : block_sample_(block_sample)
  // memory allocation
  , result_(boost::extents[block_sample->block_count()][block_sample->block_size()])
{
}

template <typename tcf_type>
void correlation<tcf_type>::compute(unsigned int level)
{
    // TODO
    // input: iterate over block_sample_->index[level]
    // output: iterate over result_[level]
}

}} // namespace observables::dynamics

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_CORRELATION_HPP */
