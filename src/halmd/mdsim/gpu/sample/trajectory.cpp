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

#include <exception>

#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/sample/trajectory.hpp>
#include <halmd/mdsim/gpu/sample/trajectory_wrapper.cuh>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace sample
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void trajectory<mdsim::samples::gpu::trajectory<dimension, float_type> >::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
    module<box_type>::required(vm);
}

template <int dimension, typename float_type>
void trajectory<mdsim::samples::host::trajectory<dimension, float_type> >::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
    module<box_type>::required(vm);
}

template <int dimension, typename float_type>
trajectory<mdsim::samples::gpu::trajectory<dimension, float_type> >::trajectory(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  , box(module<box_type>::fetch(vm))
{}

template <int dimension, typename float_type>
trajectory<mdsim::samples::host::trajectory<dimension, float_type> >::trajectory(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  , box(module<box_type>::fetch(vm))
{}

/**
 * Sample trajectory
 */
template <int dimension, typename float_type>
void trajectory<mdsim::samples::gpu::trajectory<dimension, float_type> >::acquire()
{
    // FIXME
}

/**
 * Sample trajectory
 */
template <int dimension, typename float_type>
void trajectory<mdsim::samples::host::trajectory<dimension, float_type> >::acquire()
{
    try {
        cuda::copy(particle->g_r, particle->h_r);
        cuda::copy(particle->g_image, particle->h_image);
        cuda::copy(particle->g_v, particle->h_v);
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw runtime_error("failed to copy trajectory from GPU to host");
    }

    for (size_t i = 0; i < particle->ntype; ++i) {
        r[i].reset(new position_sample_vector(particle->ntypes[i]));
        v[i].reset(new velocity_sample_vector(particle->ntypes[i]));
    }
    for (size_t i = 0; i < particle->nbox; ++i) {
        unsigned int type, tag;
        vector_type r = untagged<vector_type>(particle->h_r[i], type);
        vector_type v = untagged<vector_type>(particle->h_v[i], tag);
        vector_type image = particle->h_image[i];
        vector_type L = static_cast<vector_type>(box->length());
        // periodically extended particle position
        (*this->r[type])[tag] = r + element_prod(image, L);
        // particle velocity
        (*this->v[type])[tag] = v;
    }
}

// explicit instantiation
template class trajectory<mdsim::samples::gpu::trajectory<3, float> >;
template class trajectory<mdsim::samples::gpu::trajectory<2, float> >;
template class trajectory<mdsim::samples::host::trajectory<3, float> >;
template class trajectory<mdsim::samples::host::trajectory<2, float> >;

}}} // namespace mdsim::gpu::sample

template class module<mdsim::gpu::sample::trajectory<mdsim::samples::gpu::trajectory<3, float> > >;
template class module<mdsim::gpu::sample::trajectory<mdsim::samples::gpu::trajectory<2, float> > >;
template class module<mdsim::gpu::sample::trajectory<mdsim::samples::host::trajectory<3, float> > >;
template class module<mdsim::gpu::sample::trajectory<mdsim::samples::host::trajectory<2, float> > >;

} // namespace halmd
