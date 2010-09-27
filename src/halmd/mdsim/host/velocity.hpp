/*
 * Copyright © 2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_VELOCITY_HPP
#define HALMD_MDSIM_HOST_VELOCITY_HPP

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host
{

template <int dimension, typename float_type>
class velocity
  : public mdsim::velocity<dimension>
{
public:
    // module definitions
    typedef velocity _Self;
    typedef mdsim::velocity<dimension> _Base;
    static void options(po::options_description& desc) {}
    static void depends();
    static void select(po::options const& vm) {}

    typedef host::particle<dimension, float_type> particle_type;
    typedef typename _Base::vector_type vector_type;

    shared_ptr<particle_type> particle;

    velocity(modules::factory& factory, po::options const& vm);
    virtual ~velocity() {}
    void rescale(double factor);
    void shift(vector_type const& delta);
    void shift_rescale(vector_type const& delta, double factor);
};

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_VELOCITY_HPP */
