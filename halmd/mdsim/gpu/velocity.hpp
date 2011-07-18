/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_VELOCITY_HPP
#define HALMD_MDSIM_GPU_VELOCITY_HPP

#include <boost/make_shared.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/velocity.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
class velocity
  : public mdsim::velocity<dimension>
{
public:
    typedef mdsim::velocity<dimension> _Base;
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename _Base::vector_type vector_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    velocity(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    virtual void rescale(double factor);
    virtual void shift(vector_type const& delta);
    virtual void shift_rescale(vector_type const& delta, double factor);

private:
    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<logger_type> logger_;
    cuda::config dim_;
};

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_VELOCITY_HPP */
