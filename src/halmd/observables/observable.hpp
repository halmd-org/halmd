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

#ifndef HALMD_OBSERVABLES_OBSERVABLE_HPP
#define HALMD_OBSERVABLES_OBSERVABLE_HPP

#include <halmd/utility/options.hpp>
#include <halmd/io/statevars/writer.hpp>

namespace halmd
{

namespace observables
{

/**
 * base class for observables
 *
 * defines the interface function sample()
 */
template <int dimension>
class observable
{
public:
    // module definitions
    typedef observable _Self;
    static void options(po::options_description& desc) {};
    static void depends() {};
    static void select(po::options const& vm) {};

    typedef io::statevars::writer<dimension> writer_type;

    observable(modules::factory& factory, po::options const& vm) {};
    virtual ~observable() {}
    virtual void register_statevars(writer_type& writer) = 0;

    // sample observable and store with given time
    virtual void sample(double time) = 0;
};

} // namespace observables

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_OBSERVABLE_HPP */
