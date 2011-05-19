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

#include <boost/shared_ptr.hpp>

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
    shared_ptr<block_sample_type> block_sample
  , unsigned int factor
  , unsigned int shift
)
  // dependency injection
  : block_sample_(block_sample)
  // member initialisation
  , interval_(block_sample_->count())
{
    // setup sampling intervals for each level
    interval_.reserve(block_sample_->count());
    unsigned int i = 1;
    for (unsigned int level = 0; level < block_sample_->count(); level += 2)
    {
        interval_.push_back(i);           // even levels
        interval_.push_back(i * shift);   // odd levels
        i *= factor;
    }
}

void blocking_scheme::acquire(uint64_t step)
{
    // trigger update of input sample(s)
    on_acquire_(step);

    LOG_TRACE("[blocking_scheme] acquire current sample");

    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < block_sample_->count(); ++i) {
        if (step % interval_[i] == 0) {
            // append current sample to block
            LOG_TRACE("[blocking_scheme] append sample to block level " << i);
            block_sample_->push_back(i);
            if (block_sample_->full(i)) {
                // correlate block data with first entry
                on_correlate_block_(i);
                // discard first entry
                block_sample_->pop_front(i);
            }
        }
    }
}

void blocking_scheme::finalise(uint64_t)
{
    // iterate over all coarse-graining levels
    for (unsigned int i = 0; i < block_sample_->count(); ++i) {
        while (!block_sample_->empty(i)) {
            // correlate block data with first entry
            on_correlate_block_(i);
            // discard first entry
            block_sample_->pop_front(i);
        }
    }
}

template <typename T>
typename signal<void (uint64_t)>::slot_function_type
acquire_wrapper(shared_ptr<T> module)
{
    return bind(&T::acquire, module, _1);
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
                class_<blocking_scheme, shared_ptr<blocking_scheme> >("blocking_scheme")
                    .def(constructor<
                        shared_ptr<blocking_scheme::block_sample_type>
                      , unsigned int, unsigned int
                    >())
                    .property("acquire", &acquire_wrapper<blocking_scheme>)
                    .property("finalise", &finalise_wrapper<blocking_scheme>)
                    .def("on_acquire", &blocking_scheme::on_acquire)
                    .def("on_correlate_block", &blocking_scheme::on_correlate_block)
            ]
        ]
    ];
    return 0;
}

}} // namespace observables::dynamics

} // namespace halmd
