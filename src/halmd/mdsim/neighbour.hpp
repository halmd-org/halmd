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

#ifndef HALMD_MDSIM_NEIGHBOUR_HPP
#define HALMD_MDSIM_NEIGHBOUR_HPP

#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
class neighbour
{
public:
    // module definitions
    typedef neighbour _Self;
    static void options(po::options_description& desc) {}
    static void depends() {}
    static void select(po::options const& vm) {}

    neighbour(modules::factory& factory, po::options const& vm) {}
    virtual ~neighbour() {}
    virtual void update() = 0;
    virtual bool check() = 0;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_NEIGHBOUR_HPP */
