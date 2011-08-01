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
    shared_ptr<clock_type const> clock
  , double maximum_lag_time
  , double resolution
  , unsigned int block_size
  , unsigned int shift
  , shared_ptr<logger_type> logger
)
  // member initialisation
  : clock_(clock)
  , block_size_(block_size)
  , logger_(logger)
{
    LOG("size of coarse-graining blocks: " << block_size_);
    if (block_size_ < 2) {
        throw std::logic_error("Minimum block size is 2.");
    }

    // optionally, compute shift of shifted coarse-graining levels ('odd' levels)
    if (shift == 0) {
        shift = static_cast<unsigned int>(std::sqrt(block_size_));
    }
    assert(shift > 0);

    // report shift only if shifted blocks ('odd' levels) are enabled
    if (shift > 1) {
        LOG("coarse-graining shift: " << shift);
    }
    else {
        LOG_DEBUG("disable shifted coarse-graining blocks");
    }

    // set up sampling intervals for each level
    step_type max_interval =
        static_cast<step_type>(maximum_lag_time / resolution) / 2; // we need at least 2 data points per level
    step_type s = 1;
    while (s * shift < max_interval)
    {
        interval_.push_back(s);               // even levels
        if (shift > 1) {                      // skip if blocks would be identical
            interval_.push_back(s * shift);   // odd levels
        }
        s *= block_size_;
    }
    unsigned int block_count = interval_.size();
    LOG("number of coarse-graining blocks: " << block_count);

    // construct associated time grid
    time_.resize(boost::extents[block_count][block_size_]);
    for (unsigned int i = 0; i < block_count; ++i) {
        for (unsigned int j = 0; j < block_size_; ++j) {
            time_[i][j] = interval_[i] * j;
        }
    }
}

void blocking_scheme::sample()
{
    // trigger update of input sample(s)
    on_sample_();

    LOG_TRACE("append sample(s)");

    // check integrity of input data
    step_type step = clock_->step();
    BOOST_FOREACH(shared_ptr<block_sample_type> block_sample, block_sample_) {
        if (block_sample->timestamp() != step) {
            throw logic_error("input sample was not updated");
        }
    }

    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < interval_.size(); ++i) {
        if (step % interval_[i] == 0) {
            // append current sample to block at level 'i' for each sample type
            LOG_TRACE("append sample(s) to blocking level " << i);
            BOOST_FOREACH(shared_ptr<block_sample_type> block_sample, block_sample_) {
                block_sample->push_back(i);
            }

            // process data if block at level 'i' is full
            //
            // Checking the first blocking scheme only is sufficient,
            // since all of them are modified synchronously
            if ((*block_sample_.begin())->full(i)) {
                // call all registered correlation modules
                // and correlate block data with first entry
                LOG_TRACE("compute correlations at blocking level " << i);
                BOOST_FOREACH(shared_ptr<correlation_base> tcf, tcf_) {
                    tcf->compute(i);
                }

                // discard first entry at level 'i' for each block structure
                BOOST_FOREACH(shared_ptr<block_sample_type> block_sample, block_sample_) {
                    // check for duplicate shared_ptr to the same block sample
                    assert(block_sample->full(i));
                    block_sample->pop_front(i);
                }
            }
        }
    }
}

void blocking_scheme::finalise()
{
    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < interval_.size(); ++i) {
        // process remaining data at level 'i'
        while (!(*block_sample_.begin())->empty(i)) {
            LOG_TRACE("compute correlations at blocking level " << i);
            // call all registered correlation modules
            // and correlate block data with first entry
            BOOST_FOREACH(shared_ptr<correlation_base> tcf, tcf_) {
                tcf->compute(i);
            }

            // discard first entry at level 'i' for each block structure
            BOOST_FOREACH(shared_ptr<block_sample_type> block_sample, block_sample_) {
                block_sample->pop_front(i);
            }
        }
    }
}

blocking_scheme::slot_function_type
static wrap_sample(shared_ptr<blocking_scheme> self)
{
    return bind(&blocking_scheme::sample, self);
}

blocking_scheme::slot_function_type
static wrap_finalise(shared_ptr<blocking_scheme> self)
{
    return bind(&blocking_scheme::finalise, self);
}

static function<blocking_scheme::block_time_type const& ()>
wrap_time(shared_ptr<blocking_scheme> self)
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
                class_<blocking_scheme, shared_ptr<blocking_scheme> >("blocking_scheme_")
                    .def(constructor<
                        shared_ptr<clock_type const>
                      , double
                      , double
                      , unsigned int
                      , unsigned int
                      , shared_ptr<logger_type>
                    >())
                    .property("finalise", &wrap_finalise)
                    .property("sample", &wrap_sample)
                    .property("block_size", &blocking_scheme::block_size)
                    .property("count", &blocking_scheme::count)
                    .property("time", &wrap_time)
                    .def("add_correlation", &blocking_scheme::add_correlation)
                    .def("add_data", &blocking_scheme::add_data)
                    .def("on_sample", &blocking_scheme::on_sample)
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
