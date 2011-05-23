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

#include <halmd/io/logger.hpp>
#include <halmd/observables/dynamics/blocking_scheme.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace dynamics
{

blocking_scheme::blocking_scheme(
    unsigned int block_count
  , unsigned int block_size
  , unsigned int shift
  , double resolution
)
  // memory allocation
  : interval_(block_count)
  , time_(boost::extents[block_count][block_size])
{
    // setup sampling intervals for each level
    interval_.reserve(block_count);
    unsigned int i = 1;
    for (unsigned int level = 0; level < block_count; level += 2) {
        interval_.push_back(i);           // even levels
        interval_.push_back(i * shift);   // odd levels
        i *= block_size;
    }

    // construct associated time grid
    for (unsigned int i = 0; i < block_count; ++i) {
        for (unsigned int j = 0; j < block_size; ++j) {
            time_[i][j] = interval_[i] * j;
        }
    }
}

void blocking_scheme::sample(uint64_t step)
{
    // trigger update of input sample(s)
    on_sample_(step);

    LOG_TRACE("[blocking_scheme] append sample(s) at step " << step);

    // check integrity of input data
    BOOST_FOREACH(shared_ptr<block_sample_type> block_sample, block_sample_) {
        if (block_sample->timestamp() != step) {
            throw logic_error("input sample was not updated");
        }
//         block_sample->make_copy();
    }

    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < interval_.size(); ++i) {
        if (step % interval_[i] == 0) {
            // append current sample to block at level 'i' for each sample type
            LOG_TRACE("[blocking_scheme] append sample(s) to blocking level " << i);
            BOOST_FOREACH(shared_ptr<block_sample_type> block_sample, block_sample_) {
                block_sample->push_back(i);
            }

            // process data if block at level 'i' is full
            //
            // Checking the first blocking scheme only is sufficient,
            // since all of them are modified synchronously
            if (block_sample_.begin()->get()->full(i)) {
                // call all registered correlation modules
                // and correlate block data with first entry
                LOG_TRACE("[blocking_scheme] compute correlations at blocking level " << i);
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
}

void blocking_scheme::finalise(uint64_t step)
{
    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < interval_.size(); ++i) {
        // process remaining data at level 'i'
        while (!block_sample_.begin()->get()->empty(i)) {
            LOG_TRACE("[blocking_scheme] compute correlations at blocking level " << i);
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

template <typename T>
typename signal<void (uint64_t)>::slot_function_type
sample_wrapper(shared_ptr<T> module)
{
    return bind(&T::sample, module, _1);
}

template <typename T>
typename signal<void (uint64_t)>::slot_function_type
finalise_wrapper(shared_ptr<T> module)
{
    return bind(&T::finalise, module, _1);
}

HALMD_LUA_API int luaopen_libhalmd_observables_dynamics_blocking_scheme(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<blocking_scheme, shared_ptr<blocking_scheme> >("blocking_scheme_")
                    .def(constructor<unsigned int, unsigned int, unsigned int, double>())
                    .property("finalise", &finalise_wrapper<blocking_scheme>)
                    .property("sample", &sample_wrapper<blocking_scheme>)
                    .def("add_correlation", &blocking_scheme::add_correlation)
                    .def("add_data", &blocking_scheme::add_data)
                    .def("on_sample", &blocking_scheme::on_sample)
            ]
        ]
    ];
    return 0;
}

}} // namespace observables::dynamics

} // namespace halmd
