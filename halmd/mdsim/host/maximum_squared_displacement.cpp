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
#include <boost/bind.hpp>
#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/maximum_squared_displacement.hpp>
#include <halmd/utility/predicates/greater.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
maximum_squared_displacement<dimension, float_type>::maximum_squared_displacement(
    shared_ptr<particle_type const> particle
  , shared_ptr<box_type const> box
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  // allocate parameters
  , r0_(particle_->nparticle())
{
}

/**
 * zero maximum squared displacement
 */
template <int dimension, typename float_type>
void maximum_squared_displacement<dimension, float_type>::zero()
{
    scoped_timer_type timer(runtime_.zero);
    copy(particle_->r.begin(), particle_->r.end(), r0_.begin());
}

/**
 * compute maximum squared displacement
 */
template <int dimension, typename float_type>
float_type maximum_squared_displacement<dimension, float_type>::compute()
{
    scoped_timer_type timer(runtime_.compute);
    float_type rr_max = 0;
    for (size_t i = 0; i < particle_->nparticle(); ++i) {
        vector_type r = particle_->r[i] - r0_[i];
        box_->reduce_periodic(r);
        rr_max = max(rr_max, inner_prod(r, r));
    }
    return rr_max;
}

template <int dimension, typename float_type>
static typename signal<void ()>::slot_function_type
wrap_zero(shared_ptr<maximum_squared_displacement<dimension, float_type> > self)
{
    return bind(&maximum_squared_displacement<dimension, float_type>::zero, self);
}

template <int dimension, typename float_type>
static typename predicates::greater<float_type>::function_type
wrap_compute(shared_ptr<maximum_squared_displacement<dimension, float_type> > self)
{
    return bind(&maximum_squared_displacement<dimension, float_type>::compute, self);
}

template <int dimension, typename float_type>
void maximum_squared_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("maximum_squared_displacement_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<maximum_squared_displacement, shared_ptr<maximum_squared_displacement> >(class_name.c_str())
                    .def(constructor<
                         shared_ptr<particle_type const>
                       , shared_ptr<box_type const>
                     >())
                    .property("zero", &wrap_zero<dimension, float_type>)
                    .property("compute", &wrap_compute<dimension, float_type>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("zero", &runtime::zero)
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &maximum_squared_displacement::runtime_)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_maximum_squared_displacement(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    maximum_squared_displacement<3, double>::luaopen(L);
    maximum_squared_displacement<2, double>::luaopen(L);
#else
    maximum_squared_displacement<3, float>::luaopen(L);
    maximum_squared_displacement<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class maximum_squared_displacement<3, double>;
template class maximum_squared_displacement<2, double>;
#else
template class maximum_squared_displacement<3, float>;
template class maximum_squared_displacement<2, float>;
#endif

} // namespace host
} // namespace mdsim
} // namespace halmd
