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

#include <utility>

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/random.hpp>
#include <halmd/mdsim/velocity.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace velocity
{

template <int dimension, typename float_type>
class boltzmann
  : public mdsim::velocity<dimension>
{
public:
    typedef mdsim::velocity<dimension> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef host::random random_type;

    shared_ptr<particle_type> particle;
    shared_ptr<random_type> random;

    static void options(po::options_description& desc);
    static void resolve(po::options const& vm);
    boltzmann(po::options const& vm);
    virtual ~boltzmann() {};
    void set();

// private:
    /** assign new velocities from Gaussian distribution of width sigma,
      * return mean velocity and mean-square velocity */
    std::pair<vector_type, float_type> gaussian(float_type sigma);
    /** shift all velocities by 'v_shift' */
    void shift(vector_type const& v_shift);
    /** rescale magnitude of all velocities by factor 'scale' */
    void rescale(float_type scale);
    /** first shift, then rescale all velocities */
    void shift_rescale(vector_type const& v_shift, float_type scale);

protected:
    /** temperature */
    float_type temp_;
};

}}} // namespace mdsim::host::velocity

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_VELOCITY_BOLTZMANN_HPP */
