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

#ifndef HALMD_RNG_RANDOM_HPP
#define HALMD_RNG_RANDOM_HPP

#include <halmd/utility/options.hpp>

namespace halmd
{
namespace rng
{

class random
{
public:
    // module definitions
    typedef random _Self;
    static void options(po::options_description& desc);
    static void depends() {}
    static void select(po::options const& vm) {}

    random(po::options const& vm) {}
    virtual ~random() {}
    virtual void seed(unsigned int value) = 0;

protected:
    unsigned int readint(std::string const& fn);
    void set_seed(po::options const& vm);
};

} // namespace rng

} // namespace halmd

#endif /* ! HALMD_RNG_RANDOM_HPP */
