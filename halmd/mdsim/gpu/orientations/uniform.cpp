/*
 * Copyright © 2016       Manuel Dibak
 * Copyright © 2008-2011  Felix Höfling
 * Copyright © 2013       Nicolas Höft
 * Copyright © 2008-2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/mdsim/gpu/orientations/uniform_kernel.hpp>
#include <halmd/mdsim/gpu/orientations/uniform.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace orientations {

template <int dimension, typename float_type>
uniform<dimension, float_type>::uniform(std::shared_ptr<particle_type> particle) 
  : particle_(particle){}

template <int dimension, typename float_type>
void uniform<dimension, float_type>::set()
{
    scoped_timer_type timer(runtime_.set);
    auto orientation = make_cache_mutable(particle_->orientation());
    auto first = &*orientation->begin();
    auto last = &*orientation->end();
    size_t npart = last - first;
    try {
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        get_uniform_kernel().uniform(first, npart);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to generate uniform particle orientation on GPU");
        throw;
    }
    LOG("setting particle orientations");
}

template <int dimension, typename float_type>
void uniform<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("orientations")
            [
                class_<uniform>()
                    .def("set", &uniform::set)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &uniform::runtime_)
              , def("uniform", &std::make_shared<uniform
                  , std::shared_ptr<particle_type>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_orientations_uniform(lua_State* L)
{
    uniform<3, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class uniform<3, float>;

} // namespace mdsim
} // namespace gpu
} // namespace orientations
} // namespace halmd
