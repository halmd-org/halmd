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

#ifndef HALMD_MDSIM_HOST_VERLET_HPP
#define HALMD_MDSIM_HOST_VERLET_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/force.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/mdsim/neighbor.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim { namespace host { namespace integrator
{

template <int dimension, typename float_type>
class verlet : public mdsim::integrator<dimension, float_type>
{
public:
    typedef mdsim::integrator<dimension, float_type> _Base;
    typedef vector<float_type, dimension> vector_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef boost::shared_ptr<mdsim::particle<dimension, float_type> > particle_ptr;
    typedef mdsim::box<dimension, float_type> box_type;
    typedef boost::shared_ptr<mdsim::box<dimension, float_type> > box_ptr;
    typedef mdsim::force<dimension, float_type> force_type;
    typedef boost::shared_ptr<mdsim::force<dimension, float_type> > force_ptr;
    typedef mdsim::neighbor<dimension, float_type> neighbor_type;
    typedef boost::shared_ptr<mdsim::neighbor<dimension, float_type> > neighbor_ptr;

public:
    verlet(particle_ptr particle, box_ptr box, force_ptr force, neighbor_ptr neighbor, options const& vm);
    virtual ~verlet() {}
    void integrate(uint64_t steps);
    float_type timestep() { return timestep_; }

public:
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    boost::shared_ptr<force_type> force;
    boost::shared_ptr<neighbor_type> neighbor;

protected:
    void pre_force();
    void post_force();

protected:
    /** integration time-step */
    float_type timestep_;
    /** half time-step */
    float_type timestep_half_;
};

}}}} // namespace halmd::mdsim::host::integrator

#endif /* ! HALMD_MDSIM_HOST_VERLET_HPP */
