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

#include <boost/make_shared.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class velocity
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef fixed_vector<double, dimension> vector_type;
    typedef logger logger_type;

    velocity(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    void rescale(double factor);
    void shift(vector_type const& delta);
    void shift_rescale(vector_type const& delta, double factor);

private:
    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<logger_type> logger_;
};

} // namespace mdsim
} // namespace host
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_VELOCITY_HPP */
