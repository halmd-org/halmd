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

#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/sampler/trajectory.hpp>
#include <halmd/mdsim/gpu/sampler/trajectory_kernel.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace sampler
{

template <int dimension, typename float_type>
trajectory<mdsim::samples::gpu::trajectory<dimension, float_type> >::trajectory(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<core_type> core
)
  : _Base(particle)
  , particle(particle)  //< mdsim::host:particle
  , box(box)
  , core(core)
{
}

template <int dimension, typename float_type>
trajectory<mdsim::samples::host::trajectory<dimension, float_type> >::trajectory(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<core_type> core
)
  : _Base(particle)
  , particle(particle)
  , box(box)
  , core(core)
{
}

/**
 * Sample trajectory
 */
template <int dimension, typename float_type>
void trajectory<mdsim::samples::gpu::trajectory<dimension, float_type> >::acquire()
{
    // FIXME

    time = core->time();
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

    for (size_t i = 0; i < particle->nbox; ++i) {
        unsigned int type, tag;
        vector_type r, v;
        tie(r, type) = untagged<vector_type>(particle->h_r[i]);
        tie(v, tag) = untagged<vector_type>(particle->h_v[i]);
        vector_type image = particle->h_image[i];
        vector_type L = static_cast<vector_type>(box->length());
        // periodically extended particle position
        (*this->r[type])[tag] = r + element_prod(image, L);
        // particle velocity
        (*this->v[type])[tag] = v;
    }
    time = core->time();
}

template <typename T>
static void register_lua(lua_State* L, char const* class_name)
{
    typedef typename T::_Base _Base;
    typedef typename T::particle_type particle_type;
    typedef typename T::box_type box_type;
    typedef typename T::core_type core_type;

    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    namespace_("sampler")
                    [
                        class_<T, shared_ptr<_Base>, _Base>(class_name)
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<box_type>
                               , shared_ptr<core_type>
                            >())
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        bind(&register_lua<trajectory<mdsim::samples::gpu::trajectory<3, float> > >, _1, "trajectory_gpu_3_")
    ]
    [
        bind(&register_lua<trajectory<mdsim::samples::gpu::trajectory<2, float> > >, _1, "trajectory_gpu_2_")
    ]
    [
        bind(&register_lua<trajectory<mdsim::samples::host::trajectory<3, float> > >, _1, "trajectory_host_3_")
    ]
    [
        bind(&register_lua<trajectory<mdsim::samples::host::trajectory<2, float> > >, _1, "trajectory_host_2_")
    ];
}

}}} // namespace mdsim::gpu::sample

} // namespace halmd
