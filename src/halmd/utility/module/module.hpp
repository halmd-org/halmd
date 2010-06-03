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

#ifndef HALMD_UTILITY_MODULE_MODULE_HPP
#define HALMD_UTILITY_MODULE_MODULE_HPP

#include <exception>

#include <halmd/utility/module/exception.hpp>
#include <halmd/utility/module/demangle.hpp>
#include <halmd/utility/module/factory.hpp>
#include <halmd/utility/module/wrapper.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::logic_error;

/**
 * This template wraps static module functions.
 */
template <typename T>
class module
{
public:
    typedef rank<T> _Rank;
    typedef wrapper<T> _Wrapper;
    typedef factory::_Module_ptr _Module_ptr;
    typedef factory::_Rank_ptr _Rank_ptr;
    typedef typename _Rank::_Module_base _Base;
    typedef builder<_Base> _Base_module;
    typedef factory::_Module_map_iterator_pair _Module_map_iterator_pair;
    typedef factory::_Module_map_iterator _Module_map_iterator;

    /**
     * returns module singleton instance(s)
     */
    struct _fetch
    {
        _fetch(po::options const& vm) : vm(vm) {}
        po::options const& vm;

        /**
         * returns required or optional instance
         */
        operator shared_ptr<T>()
        {
            shared_ptr<T> result;
            _Module_map_iterator_pair range = factory::fetch(_Rank_ptr(new _Rank));

            if (range.first != range.second) {
                // ensure that range contains only a single module
                if (range.first != --range.second) {
                    throw logic_error("ambiguous dependency " + module<T>::name());
                }

                shared_ptr<_Base_module> module_ = dynamic_pointer_cast<_Base_module>(range.first->second);
                result = dynamic_pointer_cast<T>(module_->fetch(vm));
            }
            return result;
        }

        /**
         * returns many instances
         */
        operator std::vector<shared_ptr<T> >()
        {
            std::vector<shared_ptr<T> > result;
            _Module_map_iterator_pair range = factory::fetch(_Rank_ptr(new _Rank));

            for (_Module_map_iterator it = range.first; it != range.second; ++it) {
                shared_ptr<_Base_module> module_ = dynamic_pointer_cast<_Base_module>(it->second);
                result.push_back(dynamic_pointer_cast<T>(module_->fetch(vm)));
            }
            return result;
        }
    };

    static _fetch fetch(po::options const& vm)
    {
        return _fetch(vm);
    }

    /**
     * return (demangled) module name
     */
    static std::string name()
    {
        return demangled_name<T>();
    }

    /**
     * resolve required dependencies for given module
     */
    static void required(po::options const& vm)
    {
        if (!factory::resolve(_Rank_ptr(new _Rank), vm)) {
            throw unresolvable_dependency<T>("no modules available");
        }
    }

    /**
     * resolve optional dependencies for given module
     */
    static void optional(po::options const& vm)
    {
        factory::resolve(_Rank_ptr(new _Rank), vm);
    }

private:
    struct _register
    {
        _register()
        {
            factory::_register(_Rank_ptr(new _Rank), _Module_ptr(new _Wrapper));
        }
    };

    static _register register_;
};

template <typename T> typename module<T>::_register module<T>::register_;

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_MODULE_HPP */
