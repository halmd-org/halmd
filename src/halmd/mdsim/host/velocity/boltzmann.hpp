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

#ifndef HALMD_MDSIM_HOST_VELOCITY_BOLTZMANN_HPP
#define HALMD_MDSIM_HOST_VELOCITY_BOLTZMANN_HPP

#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/random.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim { namespace host { namespace velocity
{

template <int dimension, typename float_type>
class boltzmann : public mdsim::velocity<dimension, float_type>
{
public:
    typedef mdsim::velocity<dimension, float_type> _Base;
    typedef vector<float_type, dimension> vector_type;
    typedef host::particle<dimension, float_type> particle_type;
    typedef host::random random_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<random_type> random;

    boltzmann(options const& vm);
    virtual ~boltzmann() {};
    void set();

protected:
    /** temperature */
    float_type temp_;
};

}}}} // namespace halmd::mdsim::host::velocity

#endif /* ! HALMD_MDSIM_HOST_VELOCITY_BOLTZMANN_HPP */
