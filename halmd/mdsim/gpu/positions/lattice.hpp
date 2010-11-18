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

#ifndef HALMD_MDSIM_GPU_POSITIONS_LATTICE_HPP
#define HALMD_MDSIM_GPU_POSITIONS_LATTICE_HPP

#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/position.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace positions
{

template <int dimension, typename float_type, typename RandomNumberGenerator>
class lattice
  : public mdsim::position<dimension>
{
public:
    typedef mdsim::position<dimension> _Base;
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef random::gpu::random<RandomNumberGenerator> random_type;
    typedef utility::profiler profiler_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename type_traits<dimension, float>::vector_type gpu_vector_type;
    typedef typename type_traits<dimension, unsigned int>::vector_type index_type;

    static char const* module_name() { return "lattice"; }

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    boost::shared_ptr<random_type> random;

    static void luaopen(lua_State* L);

    lattice(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<random_type> random
    );
    virtual void set();
    void register_runtimes(profiler_type& profiler);

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG(set_, "setting particle positions on lattice");

private:
    boost::fusion::map<
        boost::fusion::pair<set_, accumulator<double> >
    > runtime_;
};

}}} // namespace mdsim::gpu::positions

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POSITIONS_LATTICE_HPP */
