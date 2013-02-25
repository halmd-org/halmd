/*
 * Copyright © 2011  Felix Höfling
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

#ifndef TEST_UNIT_MDSIM_POTENTIALS_GPU_NEIGHBOUR_CHAIN_HPP
#define TEST_UNIT_MDSIM_POTENTIALS_GPU_NEIGHBOUR_CHAIN_HPP

#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

/** GPU neighbour list returning a chain
 *
 * This class constructs a set of neighbour lists where each particle interacts
 * only with its successor (in the list of particles). The last particle
 * interacts with the first one.
 */

template <int dimension, typename float_type>
class neighbour_chain
  : public halmd::mdsim::gpu::neighbour
{
public:
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;

    neighbour_chain(
        std::shared_ptr<particle_type const> particle
    );

    /** neighbour lists */
    virtual halmd::cache<cuda::vector<unsigned int>> const& g_neighbour()
    {
        return g_neighbour_;
    }

    /** number of placeholders per neighbour list */
    virtual unsigned int size() const { return 1; }

    /** neighbour list stride */
    virtual unsigned int stride() const
    {
        return stride_;
    }

private:
    unsigned int stride_;
    /** neighbour lists */
    halmd::cache<cuda::vector<unsigned int>> g_neighbour_;
};

template <int dimension, typename float_type>
neighbour_chain<dimension, float_type>::neighbour_chain(
    std::shared_ptr<particle_type const> particle
)
  // member initialisation
  : stride_(particle->dim.threads())  // total number of particles and ghost particles
{
    auto g_neighbour = make_cache_mutable(g_neighbour_);

    // allocate neighbour lists of size 1
    g_neighbour->resize(stride_);

    // fill each particle's neighbour list with the following particle
    cuda::host::vector<unsigned int> neighbour(g_neighbour->size());
    for (unsigned int i = 0; i < neighbour.size(); ++i) {
        neighbour[i] = (i + 1) % particle->nparticle(); // point only to real particles
    }
    cuda::copy(neighbour.begin(), neighbour.end(), g_neighbour->begin()); // copy data to GPU
}

#endif /* ! TEST_UNIT_MDSIM_POTENTIALS_GPU_NEIGHBOUR_CHAIN_HPP */
