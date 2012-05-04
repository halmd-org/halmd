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

#include <algorithm>
#include <boost/foreach.hpp>
#include <cassert>
#include <exception>

#include <halmd/observables/dynamics/blocking_scheme.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace dynamics {

blocking_scheme::blocking_scheme(
    boost::shared_ptr<clock_type const> clock
  , double maximum_lag_time
  , double resolution
  , unsigned int block_size
  , unsigned int shift
  , unsigned int separation
  , boost::shared_ptr<logger_type> logger
)
  // member initialisation
  : clock_(clock)
  , logger_(logger)
  , block_size_(block_size)
  , separation_(separation)
  , timestep_(clock_->timestep())
{
    LOG("size of coarse-graining blocks: " << block_size_);
    if (block_size_ < 2) {
        throw invalid_argument("Minimum size of coarse-graining blocks is 2.");
    }

    // optionally, compute shift of shifted coarse-graining levels ('odd' levels)
    if (shift == 0) {
        shift = static_cast<unsigned int>(std::sqrt(block_size_));
    }
    else if (shift >= block_size_) {
        throw invalid_argument("Coarse-graining shift must be less than block size.");
    }
    assert(shift > 0);

    // report shift only if shifted blocks ('odd' levels) are enabled
    if (shift > 1) {
        LOG("coarse-graining shift: " << shift);
    }
    else {
        LOG("disable shifted coarse-graining blocks");
    }

    LOG("minimal separation of samples in time: " << separation_ * resolution);

    // set up sampling intervals for each level
    step_type max_interval = static_cast<step_type>(maximum_lag_time / clock_->timestep());
    step_type s = static_cast<step_type>(resolution / clock_->timestep());
    while (s <= max_interval)
    {
        interval_.push_back(s);               // even levels
        if (shift > 1) {                      // odd levels, skip if disabled
            step_type s_odd = s * shift;
            if (s_odd <= max_interval) {
                interval_.push_back(s_odd);
            }
        }
        s *= block_size_;
    }
    unsigned int block_count = interval_.size();
    LOG("number of coarse-graining blocks: " << block_count);
    if (block_count == 0) {
        LOG_WARNING("temporal resolution (" << resolution << ") exceeds simulation time (" << maximum_lag_time << ")");
        throw invalid_argument("Sampling interval of dynamic correlations is too large.");
    }

    // setup initial time origin for each level
    origin_.resize(block_count, 0);

    // construct associated time grid
    time_.resize(boost::extents[block_count][block_size_]);
    for (unsigned int i = 0; i < block_count; ++i) {
        for (unsigned int j = 0; j < block_size_; ++j) {
            time_[i][j] = (interval_[i] * j) * clock_->timestep();
        }
    }
}

connection blocking_scheme::on_correlate(boost::shared_ptr<correlation_base> tcf)
{
    assert(find(tcf_.begin(), tcf_.end(), tcf) == tcf_.end());
    return tcf_.connect(tcf);
}

connection blocking_scheme::on_sample(boost::shared_ptr<block_sample_type> block_sample)
{
    assert(block_sample->count() == count());
    assert(block_sample->block_size() == block_size_);
    // check for duplicate shared_ptr to the same block sample,
    // which would result in duplicate push and pop operations
    assert(find(block_sample_.begin(), block_sample_.end(), block_sample) == block_sample_.end());
    return block_sample_.connect(block_sample);
}

void blocking_scheme::sample()
{
    // FIXME only if needed
    // trigger update of input sample(s)
    on_prepend_sample_();

    LOG_TRACE("append sample(s)");

    // check time step is same as upon construction
    if (clock_->timestep() != timestep_) {
        throw logic_error("blocking scheme does not allow variable time step");
    }

    // check integrity of input data
    step_type step = clock_->step();
    BOOST_FOREACH(boost::shared_ptr<block_sample_type> block_sample, block_sample_) {
        if (block_sample->step() != step) {
            throw logic_error("input sample was not updated");
        }
    }

    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < interval_.size(); ++i) {
        if (step % interval_[i] == 0 && step >= origin_[i]) {
            // append current sample to block at level 'i' for each sample type
            LOG_TRACE("append sample(s) to blocking level " << i);
            BOOST_FOREACH(boost::shared_ptr<block_sample_type> block_sample, block_sample_) {
                block_sample->push_back(i);
            }

            // process data if block at level 'i' is full
            //
            // Checking the first blocking scheme only is sufficient,
            // since all of them are modified synchronously
            if (!block_sample_.empty() && (*block_sample_.begin())->full(i)) {
                process(i);
            }
        }
    }
    on_append_sample_();
}

void blocking_scheme::finalise()
{
    on_prepend_finalise_();
    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < interval_.size(); ++i) {
        // process remaining data at level 'i'
        while (!block_sample_.empty() && !(*block_sample_.begin())->empty(i)) {
            process(i);
        }
    }
    on_append_finalise_();
}

void blocking_scheme::process(unsigned int level)
{
    // call all registered correlation modules
    // and correlate block data with first entry
    LOG_TRACE("compute correlations at blocking level " << level << " from step " << origin_[level]);
    BOOST_FOREACH(boost::shared_ptr<correlation_base> tcf, tcf_) {
        tcf->compute(level);
    }

    // update time origin for next computation at this level
    //
    // make sure that the new origin is a multiple of this level's sampling interval
    // and at least incremented by one interval
    unsigned int skip = max((separation_ + interval_[level] - 1) / interval_[level], step_type(1));
    origin_[level] += skip * interval_[level];

    // for each block structure, discard all entries earlier
    // than the new time origin at this level
    LOG_TRACE("discard first " << skip << " entries at level " << level);
    BOOST_FOREACH(boost::shared_ptr<block_sample_type> block_sample, block_sample_) {
        for (unsigned int k = 0; k < skip && !block_sample->empty(level); ++k) {
            block_sample->pop_front(level);
        }
    }
}

connection blocking_scheme::on_prepend_sample(slot_function_type const& slot)
{
    return on_prepend_sample_.connect(slot);
}

connection blocking_scheme::on_append_sample(slot_function_type const& slot)
{
    return on_append_sample_.connect(slot);
}

connection blocking_scheme::on_prepend_finalise(slot_function_type const& slot)
{
    return on_prepend_finalise_.connect(slot);
}

connection blocking_scheme::on_append_finalise(slot_function_type const& slot)
{
    return on_append_finalise_.connect(slot);
}

blocking_scheme::slot_function_type
static wrap_sample(boost::shared_ptr<blocking_scheme> self)
{
    return bind(&blocking_scheme::sample, self);
}

blocking_scheme::slot_function_type
static wrap_finalise(boost::shared_ptr<blocking_scheme> self)
{
    return bind(&blocking_scheme::finalise, self);
}

static boost::function<blocking_scheme::block_time_type const& ()>
wrap_time(boost::shared_ptr<blocking_scheme> self)
{
    return bind(&blocking_scheme::time, self);
}

void blocking_scheme::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<blocking_scheme, boost::shared_ptr<blocking_scheme> >("blocking_scheme")
                    .def(constructor<
                        boost::shared_ptr<clock_type const>
                      , double
                      , double
                      , unsigned int
                      , unsigned int
                      , unsigned int
                      , boost::shared_ptr<logger_type>
                    >())
                    .property("finalise", &wrap_finalise)
                    .property("sample", &wrap_sample)
                    .property("block_size", &blocking_scheme::block_size)
                    .property("separation", &blocking_scheme::separation)
                    .property("count", &blocking_scheme::count)
                    .property("time", &wrap_time)
                    .def("on_correlate", &blocking_scheme::on_correlate)
                    .def("on_sample", &blocking_scheme::on_sample)
                    .def("on_prepend_sample", &blocking_scheme::on_prepend_sample)
                    .def("on_append_sample", &blocking_scheme::on_append_sample)
                    .def("on_prepend_finalise", &blocking_scheme::on_prepend_finalise)
                    .def("on_append_finalise", &blocking_scheme::on_append_finalise)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_dynamics_blocking_scheme(lua_State* L)
{
    blocking_scheme::luaopen(L);
    return 0;
}

} // namespace dynamics
} // namespace observables
} // namespace halmd
