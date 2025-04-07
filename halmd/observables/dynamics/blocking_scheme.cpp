/*
 * Copyright © 2011  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>

#include <halmd/observables/dynamics/blocking_scheme.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace observables {
namespace dynamics {

blocking_scheme::blocking_scheme(
    std::shared_ptr<clock_type const> clock
  , time_type maximum_lag_time
  , time_type resolution
  , unsigned int block_size
  , unsigned int shift
  , unsigned int separation
  , std::shared_ptr<logger> logger
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
    if (shift == 0) {
        throw invalid_argument("Coarse-graining shift must be non-zero.");
    }
    if (shift >= block_size_) {
        throw invalid_argument("Coarse-graining shift must be less than block size.");
    }

    // report shift only if shifted blocks ('odd' levels) are enabled
    if (shift > 1) {
        LOG("coarse-graining shift: " << shift);
    }
    else {
        LOG("disable shifted coarse-graining blocks");
    }

    // convert resolution and max lag time to steps using rounding to avoid truncation, e.g.,
    // avoid that s = 0 if resolution ≈ timestep_ up to floating-point precision
    if (resolution <= 0) {
        throw invalid_argument("Resolution must be a positive value.");
    }
    step_type max_interval = static_cast<step_type>(std::round(maximum_lag_time / timestep_));
    step_type s = static_cast<step_type>(std::round(resolution / timestep_));
    resolution = s * timestep_;

    LOG("minimal separation of samples in time: " << separation_ * resolution);

    // set up sampling intervals for each level
    assert(s > 0);
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
    origin_.resize(block_count, clock_->step());

    // construct associated time grid
    time_.resize(boost::extents[block_count][block_size_]);
    for (unsigned int i = 0; i < block_count; ++i) {
        for (unsigned int j = 0; j < block_size_; ++j) {
            time_[i][j] = (interval_[i] * j) * clock_->timestep();
        }
    }
}

connection blocking_scheme::on_correlate(std::shared_ptr<correlation_base> tcf)
{
    assert(find(tcf_.begin(), tcf_.end(), tcf) == tcf_.end());
    return tcf_.connect(tcf);
}

connection blocking_scheme::on_sample(std::shared_ptr<block_sample_type> block_sample)
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
    on_prepend_sample_();

    LOG_DEBUG("append sample(s)");

    // check time step is same as upon construction
    if (clock_->timestep() != timestep_) {
        throw logic_error("blocking scheme does not allow variable time step");
    }

    step_type const step = clock_->step();

    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < interval_.size(); ++i) {
        if ((step - origin_[i]) % interval_[i] == 0 && step >= origin_[i]) {
            // append current sample to block at level 'i' for each sample type
            LOG_TRACE("append sample(s) to blocking level " << i);
            for (std::shared_ptr<block_sample_type> block_sample : block_sample_) {
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
    LOG_DEBUG("compute correlations at blocking level " << level << " from step " << origin_[level]);
    for (std::shared_ptr<correlation_base> tcf : tcf_) {
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
    for (std::shared_ptr<block_sample_type> block_sample : block_sample_) {
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

static std::function<void ()>
wrap_sample(std::shared_ptr<blocking_scheme> self)
{
    return [=]() {
        self->sample();
    };
}

static std::function<void ()>
wrap_finalise(std::shared_ptr<blocking_scheme> self)
{
    return [=]() {
        self->finalise();
    };
}

static std::function<blocking_scheme::block_time_type const& ()>
wrap_time(std::shared_ptr<blocking_scheme> self)
{
    return [=]() -> blocking_scheme::block_time_type const& {
        return self->time();
    };
}

void blocking_scheme::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<blocking_scheme, std::shared_ptr<blocking_scheme> >("blocking_scheme")
                    .def(constructor<
                        std::shared_ptr<clock_type const>
                      , time_type
                      , time_type
                      , unsigned int
                      , unsigned int
                      , unsigned int
                      , std::shared_ptr<logger>
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
