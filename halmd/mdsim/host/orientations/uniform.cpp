/*
 * Copyright © 2017       Jake Atwell
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

#include <halmd/mdsim/host/orientations/uniform.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace orientations {

template <int dimension, typename float_type>
uniform<dimension, float_type>::uniform(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<random_type> random
  , std::shared_ptr<logger> logger
) 
  : particle_(particle)
  , random_(random)
  , logger_(logger)
{
}

template <int dimension, typename float_type>
void uniform<dimension, float_type>::set()
{
    scoped_timer_type timer(runtime_.set);
    auto orientation = make_cache_mutable(particle_->orientation());

    //assigns a random unit-vector efficiently for each particle
    for (auto& u : *orientation)
    {
        random_->unit_vector(u); 
        
    }
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
                  , std::shared_ptr<random_type>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_orientations_uniform(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    //uniform<2, double>::luaopen(L);
    uniform<3, double>::luaopen(L);
#else    
    //uniform<2, float>::luaopen(L);
    uniform<3, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class uniform<2, double>;
template class uniform<3, double>;
#else
template class uniform<2, float>;
template class uniform<3, float>;
#endif

} // namespace mdsim
} // namespace host 
} // namespace orientations
} // namespace halmd
