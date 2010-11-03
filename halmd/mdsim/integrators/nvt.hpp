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

#ifndef HALMD_MDSIM_INTEGRATORS_NVT_HPP
#define HALMD_MDSIM_INTEGRATORS_NVT_HPP

#include <lua.hpp>

#include <halmd/mdsim/integrator.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim { namespace integrators
{

template <int dimension>
class nvt : public integrator<dimension>
{
public:
    typedef mdsim::integrator<dimension> _Base;

    static void options(po::options_description& desc);
    static void luaopen(lua_State* L);

    nvt() {}
    virtual ~nvt() {}
    virtual double temperature() const = 0;
    virtual void temperature(double temperature) = 0;
};

}} // namespace mdsim::integrators

} // namespace halmd

#endif /* ! HALMD_MDSIM_INTEGRATORS_NVT_HPP */
