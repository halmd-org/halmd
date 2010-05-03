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

#include <boost/logic/tribool.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <set>
#include <typeinfo>

#include <halmd/utility/module/builder.hpp>
#include <halmd/utility/module/exception.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace utility { namespace module
{

// import into namespace
using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using boost::weak_ptr;

/**
 * Concrete module
 */
template <typename T>
class module
  : public builder<T>
{
public:
    typedef typename builder<T>::_Module_base _Base;

    module() : resolved_(boost::indeterminate) {}

    /**
     * weak module ordering
     */
    bool rank(shared_ptr<builder<> > const& other) const
    {
        // For the case that the *other* module derives from
        // *this* module and should thus be ranked higher than
        // this module, returns false. Otherwise returns true.
        return !dynamic_pointer_cast<builder<T> >(other);
    }

    /**
     * returns singleton instance
     */
    shared_ptr<_Base> fetch(po::options const& vm)
    {
        // We use an observing weak pointer instead of an owning
        // shared pointer to let the caller decide when the
        // singleton instance and its dependencies are destroyed.
        //
        // Special care has to be taken not to destroy the
        // instance before returning it over to the caller.

        shared_ptr<T> singleton(singleton_.lock());
        if (!singleton) {
            singleton.reset(new T(vm));
            singleton_ = singleton;
        }
        return singleton;
    }

    /**
     * assemble module options
     */
    void options(po::options_description& desc)
    {
        builder<T>::options(desc);
    }

    /**
     * resolve module dependencies
     */
    void resolve(po::options const& vm)
    {
        if (resolved_) {
            LOG_DEBUG("cached resolvable module");
        }
        else if (!resolved_) {
            LOG_DEBUG("cached irresolvable module");
            throw module_exception("irresolvable module " + name());
        }
        else {
            resolved_ = false;
            builder<T>::resolve(vm);
            resolved_ = true;
        }
    }

    /**
     * returns module name
     */
    std::string name()
    {
        return typeid(T).name();
    }

    /** module instance observer */
    static weak_ptr<T> singleton_;

private:
    boost::tribool resolved_;
};

template <typename T> weak_ptr<T> module<T>::singleton_;

}} // namespace utility::module

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_MODULE_HPP */
