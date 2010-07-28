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

#ifndef HALMD_MDSIM_FORCE_HPP
#define HALMD_MDSIM_FORCE_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/mdsim/particle.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim
{

/**
 * The force module computes all interparticle forces.
 * The current particle positions are read from and
 * the result is stored in the particle module. Periodic
 * boundaries may be taken into account by reference to the
 * box module.
 */

template <int dimension>
class force
{
public:
    // module definitions
    typedef force _Self;
    static void options(po::options_description& desc);
    static void depends();
    static void select(po::options const& vm) {}

    typedef mdsim::particle<dimension> particle_type;

    shared_ptr<particle_type> particle;

    force(modules::factory& factory, po::options const& vm);
    virtual ~force() {}
    virtual void compute() = 0;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_FORCE_HPP */
