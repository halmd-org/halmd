/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <boost/bind.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/module/exception.hpp>
#include <halmd/utility/module/factory.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace utility { namespace module
{

/**
 * register module builder at program startup
 */
void factory::_register(_Rank_ptr rank_, _Module_ptr module_)
{
    if (!modules().insert(make_pair(rank_, module_)).second) {
        throw module_error("duplicate module " + module_->name());
    }
}

/**
 * resolve module dependencies
 *
 * returns the number of resolved modules, which is used to
 * validate a dependency as required, optional or one-to-many.
 */
size_t factory::resolve(_Rank_ptr rank_, po::options const& vm)
{
    // Given a required, optional or one-to-many dependency,
    // which may be a (multiply) derived or base class, this
    // function looks for suitable modules.
    // A module is suitable (1) if it is or derives from the
    // given dependency type and (2) if its dependencies can be
    // resolved.
    //

#define STACK_DEBUG(x) LOG_DEBUG("[" << stack_.size() << "] " << x)

    // We cache the result of a dependency resolution by holding
    // the number of resolved modules in a rank-indexed map.

    if (!cache_.count(rank_)) {
        STACK_DEBUG("resolve dependency " << rank_->name());

        set<_Rank_ptr, rank_order_equal_base> resolved;

        // Check each of the modules registered in the base class
        // factory during program startup for suitability.

        _Module_map_iterator_pair range = fetch(rank_);

        for (_Module_map_iterator it = range.first; it != range.second; ) {
            // Check if the set of suitable modules contains a
            // module which derives from this module. The module
            // set ordering guarantees that derived modules
            // come before base modules in the sequence.

            if (resolved.count(it->first)) {
                STACK_DEBUG("skipping base module " << it->second->name());
                modules().erase(it++);
                continue;
            }

            STACK_DEBUG("resolve module " << it->second->name());

            // Take note of the current top of the dependency
            // resolution stack to rewind in case of failure.
            _Module_stack_iterator top = stack_.end();
            stack_.push_back(it->second);
            // Try to resolve the dependencies of the module.
            try {
                it->second->resolve(vm);
            }
            catch (module_error const& e) {
                // The module is irresolvable, therefore we
                // rewind the stack and erase this module.
                STACK_DEBUG(e.what());
                stack_.erase(top, stack_.end());
                modules().erase(it++);
                continue;
            }
            // The module is resolvable.
            resolved.insert(it->first);
            ++it;
        }
        cache_.insert(make_pair(rank_, resolved.size()));
    }

#undef STACK_DEBUG

    return cache_.at(rank_);
}

/**
 * assemble module options
 */
po::options_description factory::options()
{
    po::options_description desc;

    // Using the stack of modules created during dependency
    // resolution we build a unique set of used modules and
    // collect the options of each module.

    set<_Module_ptr> set(stack_.begin(), stack_.end());
    for_each(set.begin(), set.end(), bind(&builder<>::options, _1, ref(desc)));
    return desc;
}

/**
 * module map equality for modules with equal or derived dank
 */
struct derived_rank_equality
{
    /** search rank */
    typedef factory::_Rank_ptr _Rank_ptr;
    /** module rank */
    typedef factory::_Module_map::value_type _Module_pair;

    bool operator()(_Module_pair left, _Rank_ptr right) const
    {
        // If the search rank is right in the inequality and the
        // module rank is left, treat equal or derived module
        // rank as equal to the search rank.

        return rank_order_equal_base()(left.first, right);
    }

    bool operator()(_Rank_ptr left, _Module_pair right) const
    {
        // If the search rank is left in the inequality and the
        // module rank is right, use normal rank ordering. We are
        // looking only for equal or derived module ranks, which
        // are ordered left to base module ranks in the map.

        return rank_order()(left, right.first);
    }
};

/**
 * returns a range of modules with equal or derived rank
 */
factory::_Module_map_iterator_pair factory::fetch(_Rank_ptr rank_)
{
    return equal_range(modules().begin(), modules().end(), rank_, derived_rank_equality());
}

/**
 * returns singleton builder set
 */
factory::_Module_map& factory::modules()
{
    // What's the "static initialization order fiasco"?
    // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12

    static _Module_map modules_;
    return modules_;
}

/** stack to keep track of used modules */
factory::_Module_stack factory::stack_;
/** resolved module cache */
factory::_Rank_cache factory::cache_;

}} // namespace utility::module

} // namespace halmd
