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

#ifndef HALMD_UTILITY_MODULE_FACTORY_HPP
#define HALMD_UTILITY_MODULE_FACTORY_HPP

#include <boost/shared_ptr.hpp>
#include <deque>
#include <map>
#include <set>

#include <halmd/utility/module/builder.hpp>
#include <halmd/utility/module/rank.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::shared_ptr;

/**
 * The factory manages available modules, which are registered
 * automagically during program startup. To provide a module with
 * singleton instantiations of its dependent modules, the factory
 * has a sophisticated dependency resolution algorithm. This
 * algorithm is also used to keep track of used modules for
 * assembling program options.
 */
class factory
{
public:
    typedef shared_ptr<rank<> > _Rank_ptr;
    typedef shared_ptr<builder<> > _Module_ptr;
    typedef std::map<_Rank_ptr, _Module_ptr, rank_order> _Module_map;
    typedef _Module_map::iterator _Module_map_iterator;
    typedef std::pair<_Module_map_iterator, _Module_map_iterator> _Module_map_iterator_pair;
    typedef std::map<_Rank_ptr, size_t, rank_order> _Rank_cache;
    typedef std::deque<_Module_ptr> _Module_stack;
    typedef _Module_stack::iterator _Module_stack_iterator;

    static void _register(_Rank_ptr rank_, _Module_ptr module_);
    static _Module_map_iterator_pair fetch(_Rank_ptr rank_);
    static size_t resolve(_Rank_ptr rank_, po::options const& vm);
    static po::options_description options();

private:
    static _Module_map& modules();
    static _Rank_cache cache_;
    static _Module_stack stack_;
};

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_FACTORY_HPP */
