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

#ifndef HALMD_MDSIM_HOST_VELOCITY_FILE_HPP
#define HALMD_MDSIM_HOST_VELOCITY_FILE_HPP

#include <vector>

#include <halmd/io/trajectory/reader.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace velocity
{

template <int dimension, typename float_type>
class file
  : public mdsim::velocity<dimension>
{
public:
    // module definitions
    typedef file _Self;
    typedef mdsim::velocity<dimension> _Base;
    static void options(po::options_description& desc) {}
    static void depends();
    static void select(po::options const& vm);

    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef io::trajectory::reader<dimension> reader_type;
    typedef samples::host::trajectory<dimension, float_type> sample_type;

    shared_ptr<reader_type> reader;
    shared_ptr<sample_type> sample;
    shared_ptr<particle_type> particle;

    file(modules::factory& factory, po::options const& vm);
    virtual ~file() {}
    void set();
};

}}} // namespace mdsim::host::velocity

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_VELOCITY_FILE_HPP */
