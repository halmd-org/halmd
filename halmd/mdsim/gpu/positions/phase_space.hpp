/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_POSITIONS_PHASE_SPACE_HPP
#define HALMD_MDSIM_GPU_POSITIONS_PHASE_SPACE_HPP

#include <boost/make_shared.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {

template <int dimension, typename float_type>
class phase_space
  : public mdsim::position<dimension>
{
public:
    typedef mdsim::position<dimension> _Base;
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef observables::host::samples::phase_space<dimension, float_type> sample_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    phase_space(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<sample_type const> sample
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    virtual void set();

private:
    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<sample_type const> sample_;
    boost::shared_ptr<logger_type> logger_;
};

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POSITIONS_PHASE_SPACE_HPP */
