/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_OBSERVABLES_GPU_SAMPLES_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_GPU_SAMPLES_PHASE_SPACE_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <limits>
#include <lua.hpp>
#include <stdint.h>
#include <vector>

#include <halmd/mdsim/type_traits.hpp>

namespace halmd
{
namespace observables { namespace gpu { namespace samples
{

template <int dimension, typename float_type>
class phase_space
{
public:
    typedef typename mdsim::type_traits<dimension, float_type>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float_type>::gpu::coalesced_vector_type gpu_vector_type;

    /** sample vector type for all particles of a species */
    typedef cuda::vector<gpu_vector_type> sample_vector;
    /** sample pointer type for all particle of a species */
    typedef boost::shared_ptr<sample_vector> sample_vector_ptr;
    /** sample pointer type for all species */
    typedef std::vector<sample_vector_ptr> sample_vector_ptr_vector;

    /** periodically extended particle positions */
    sample_vector_ptr_vector r;
    /** particle velocities */
    sample_vector_ptr_vector v;
    /** simulation step when sample was taken */
    uint64_t step;

    static void luaopen(lua_State* L);

    phase_space(std::vector<unsigned int> ntypes);

    /** free shared pointers and re-allocate memory, values are not initialised */
    void reset();
};

template <int dimension, typename float_type>
inline phase_space<dimension, float_type>::phase_space(std::vector<unsigned int> ntypes)
  // allocate sample pointers
  : r(ntypes.size())
  , v(ntypes.size())
  // initialise attributes
  , step(std::numeric_limits<uint64_t>::max())
{
    for (size_t i = 0; i < ntypes.size(); ++i) {
        r[i].reset(new sample_vector(ntypes[i]));
        v[i].reset(new sample_vector(ntypes[i]));
    }
}

template <int dimension, typename float_type>
inline void phase_space<dimension, float_type>::reset()
{
    // free shared pointers and re-allocate memory
    for (size_t i = 0; i < r.size(); ++i) {
        r[i].reset(new sample_vector(r[i]->size()));
    }
    for (size_t i = 0; i < v.size(); ++i) {
        v[i].reset(new sample_vector(v[i]->size()));
    }
    // make time stamp invalid
    step = std::numeric_limits<uint64_t>::max();
}

}}} // namespace observables::gpu::samples

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_SAMPLES_PHASE_SPACE_HPP */
