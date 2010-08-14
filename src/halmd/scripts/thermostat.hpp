/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_THERMOSTAT_HPP
#define HALMD_THERMOSTAT_HPP

#include <halmd/script.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/module.hpp>

namespace halmd
{
namespace scripts
{

/**
 * Pseudo-thermostat run
 *
 * This script uses the boltzmann module to periodically assign
 * Boltzmann-distributed velocities to all particles, which is
 * sufficient to realise pre-equilibration cooling.
 */
template <int dimension>
class thermostat
  : public script<dimension>
{
public:
    // module definitions
    typedef thermostat _Self;
    typedef halmd::script<dimension> _Base;
    static void depends();
    static void options(po::options_description& desc);
    static void select(po::options const& vm) {}

    typedef typename _Base::profile_writer_type profile_writer_type;

    thermostat(modules::factory& factory, po::options const& vm);
    void run();

    using _Base::core;
    using _Base::profile_writers;

protected:
    /** heat bath collision rate */
    float rate_;
    /** heat bath coupling interval */
    unsigned int interval_;
};

} // namespace scripts

} // namespace halmd

#endif /* ! HALMD_THERMOSTAT_HPP */
