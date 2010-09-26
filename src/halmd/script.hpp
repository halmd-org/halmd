/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_SCRIPT_HPP
#define HALMD_SCRIPT_HPP

#include <halmd/utility/lua/lua_include.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{

/**
 * NVE ensemble run
 */
template <int dimension>
class script
{
public:
    script(modules::factory& factory, po::options const& vm);
    virtual ~script() {}
    virtual void load_wrapper();
    virtual void load_library();
    virtual void run();

protected:
    /** Lua state */
    boost::shared_ptr<lua_State> L_;
};

} // namespace halmd

#endif /* ! HALMD_SCRIPT_HPP */
