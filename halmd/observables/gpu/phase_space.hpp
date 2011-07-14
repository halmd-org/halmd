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

#ifndef HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP

#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/observables/phase_space.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <typename sample_type>
class phase_space;

/**
 * Sample phase_space from GPU memory to GPU memory
 */
template <int dimension, typename float_type>
class phase_space<gpu::samples::phase_space<dimension, float_type> >
  : public observables::phase_space<dimension>
{
public:
    typedef observables::phase_space<dimension> _Base;
    typedef gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;

    boost::shared_ptr<sample_type> sample;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void luaopen(lua_State* L);

    phase_space(
        boost::shared_ptr<sample_type> sample
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    virtual void acquire(uint64_t step);
};

/**
 * Sample phase_space from GPU memory to host memory
 */
template <int dimension, typename float_type>
class phase_space<host::samples::phase_space<dimension, float_type> >
  : public observables::phase_space<dimension>
{
public:
    typedef observables::phase_space<dimension> _Base;
    typedef host::samples::phase_space<dimension, float_type> sample_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef fixed_vector<float_type, dimension> vector_type;

    boost::shared_ptr<sample_type> sample;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void luaopen(lua_State* L);

    phase_space(
        boost::shared_ptr<sample_type> sample
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    virtual void acquire(uint64_t step);
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP */
