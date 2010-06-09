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

#ifndef HALMD_MDSIM_HOST_FORCE_HPP
#define HALMD_MDSIM_HOST_FORCE_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/host/box.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/mdsim/host/neighbour.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/thermodynamics.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host
{

template <int dimension, typename float_type>
class force
  : public mdsim::force<dimension>
{
public:
    typedef mdsim::force<dimension> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::matrix_type matrix_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef host::box<dimension> box_type;
    typedef host::neighbour<dimension, float_type> neighbour_type;
    typedef host::thermodynamics<dimension> thermodynamics_type;
    typedef host::forces::smooth<dimension, float_type> smooth_type;

    shared_ptr<particle_type> particle;
    shared_ptr<box_type> box;
    shared_ptr<neighbour_type> neighbour;
    shared_ptr<thermodynamics_type> thermodynamics;
    shared_ptr<smooth_type> smooth;

    static void options(po::options_description& desc);
    static void resolve(po::options const& vm);
    force(po::options const& vm);
    virtual ~force() {};
    virtual void compute() = 0;
    virtual matrix_type const& cutoff() = 0;
};

}} // namespace mdsim::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCE_HPP */
