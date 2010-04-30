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

#ifndef HALMD_MDSIM_HOST_SAMPLE_HPP
#define HALMD_MDSIM_HOST_SAMPLE_HPP

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/numeric/host/blas/vector.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace sample
{

template <int dimension, typename float_type>
class trajectory
  : public mdsim::samples::host::trajectory<dimension, float_type>
{
public:
    typedef mdsim::samples::host::trajectory<dimension, float_type> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void options(po::options_description& desc) {}
    static void resolve(po::options const& vm);
    trajectory(po::options const& vm);
    virtual ~trajectory() {}
    void acquire();

    typedef typename _Base::position_sample_vector position_sample_vector;
    typedef typename _Base::velocity_sample_vector velocity_sample_vector;

    /** periodically extended particle positions */
    using _Base::r;
    /** particle velocities */
    using _Base::v;
};

}}} // namespace mdsim::host::sample

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_SAMPLE_HPP */
