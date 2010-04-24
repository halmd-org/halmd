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

#ifndef HALMD_MDSIM_POSITION_HPP
#define HALMD_MDSIM_POSITION_HPP

#include <halmd/utility/module.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
class position
{
public:
    typedef typename boost::shared_ptr<position> module_ptr;

    position(options const& vm) {}
    virtual ~position() {}
    virtual void set() = 0;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_POSITION_HPP */
