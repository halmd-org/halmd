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

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/velocities/phase_space.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    shared_ptr<particle_type> particle
  , shared_ptr<sample_type const> sample
  , shared_ptr<logger_type> logger
)
  : _Base(particle, logger)
  // dependency injection
  , particle_(particle)
  , sample_(sample)
  , logger_(logger)
{
}

/**
 * set particle velocities
 */
template <int dimension, typename float_type>
void phase_space<dimension, float_type>::set()
{
    LOG("set particle velocities from phase space sample");

    scoped_timer_type timer(runtime_.set);

    // assign particle velocities and tags
    size_t n = 0; // indicates the boundary to the next particle type
    for (size_t j = 0, i = 0; j < particle_->ntype; ++j) {
        typename sample_type::sample_vector const& v_sample = *sample_->v[j];
        n += particle_->ntypes[j];
        assert(particle_->ntypes[j] == v_sample.size());
        assert(n <= particle_->h_v.size());
        for (size_t k = 0; i < n; ++i, ++k) {
            particle_->h_v[i] = particle_kernel::tagged<vector_type>(v_sample[k], k);
        }
    }

    try {
#ifdef USE_VERLET_DSFUN
        // erase particle velocity vectors (double-single precision)
        cuda::memset(particle_->g_v, 0, particle_->g_v.capacity());
#endif
        cuda::copy(particle_->h_v, particle_->g_v);
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("failed to copy particle velocities to GPU");
        throw;
    }
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("velocities")
                [
                    class_<phase_space, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                             shared_ptr<particle_type>
                           , shared_ptr<sample_type const>
                           , shared_ptr<logger_type>
                        >())
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("set", &runtime::set)
                        ]
                        .def_readonly("runtime", &phase_space::runtime_)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_velocities_phase_space(lua_State* L)
{
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class phase_space<3, float>;
template class phase_space<2, float>;

} // namespace mdsim
} // namespace gpu
} // namespace velocities
} // namespace halmd
