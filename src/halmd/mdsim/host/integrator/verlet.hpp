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
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/integrator.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim { namespace host { namespace integrator
{

template <int dimension, typename float_type>
class verlet : public mdsim::integrator<dimension>
{
public:
    typedef mdsim::integrator<dimension> _Base;
    typedef vector<float_type, dimension> vector_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;

    verlet(options const& vm);
    virtual ~verlet() {}
    void integrate();
    void finalize();
    double timestep() { return timestep_; }

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

protected:
    /** integration time-step */
    double timestep_;
    /** half time-step */
    double timestep_half_;
};

}}}} // namespace halmd::mdsim::host::integrator

#endif /* ! HALMD_MDSIM_HOST_VERLET_HPP */
