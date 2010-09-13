/*
 * Copyright © 2010  Felix Höfling
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

#ifndef HALMD_MDSIM_SAMPLER_HPP
#define HALMD_MDSIM_SAMPLER_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/io/trajectory/writer.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
class sampler
{
public:
    // module definitions
    typedef sampler _Self;
    static void options(po::options_description& desc);
    static void depends();
    static void select(po::options const& vm) {};

    typedef mdsim::core<dimension> core_type;
    typedef observables::thermodynamics<dimension> thermodynamics_type;
    typedef io::trajectory::writer<dimension> trajectory_writer_type;

    sampler(modules::factory& factory, po::options const& vm);
    void sample(bool force=false);

    shared_ptr<core_type> core;
    shared_ptr<thermodynamics_type> thermodynamics;
    shared_ptr<trajectory_writer_type> trajectory_writer;

private:
    // value from option --sampling-stat-vars
    unsigned stat_vars_interval_;
    // value from option --sampling-trajectory
    unsigned trajectory_interval_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_SAMPLER_HPP */
