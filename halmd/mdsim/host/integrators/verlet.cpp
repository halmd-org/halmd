/*
 * Copyright © 2008-2012  Peter Colberg
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
#include <cmath>
#include <memory>

#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , double timestep
  , std::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
{
    set_timestep(timestep);
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::set_timestep(double timestep)
{
    timestep_ = timestep;
    timestep_half_ = 0.5 * timestep;
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    cache_proxy<position_array_type> position = particle_->position();
    cache_proxy<image_array_type> image = particle_->image();
    cache_proxy<velocity_array_type> velocity = particle_->velocity();
    cache_proxy<force_array_type const> force = particle_->force();
    cache_proxy<mass_array_type const> mass = particle_->mass();
    size_type const nparticle = particle_->nparticle();

    scoped_timer_type timer(runtime_.integrate);

    for (size_type i = 0; i < nparticle; ++i) {
        vector_type& v = (*velocity)[i];
        vector_type& r = (*position)[i];
        v += (*force)[i] * timestep_half_ / (*mass)[i];
        r += v * timestep_;
        (*image)[i] += box_->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::finalize()
{
    cache_proxy<velocity_array_type> velocity = particle_->velocity();
    cache_proxy<force_array_type const> force = particle_->force();
    cache_proxy<mass_array_type const> mass = particle_->mass();
    size_type const nparticle = particle_->nparticle();

    scoped_timer_type timer(runtime_.finalize);

    for (size_type i = 0; i < nparticle; ++i) {
        (*velocity)[i] += (*force)[i] * timestep_half_ / (*mass)[i];
    }
}

template <typename integrator_type>
static std::function<void ()>
wrap_integrate(std::shared_ptr<integrator_type> self)
{
    return [=]() {
        self->integrate();
    };
}

template <typename integrator_type>
static std::function<void ()>
wrap_finalize(std::shared_ptr<integrator_type> self)
{
    return [=]() {
        self->finalize();
    };
}

template <typename integrator_type>
static std::function<void (double)>
wrap_set_timestep(std::shared_ptr<integrator_type> self)
{
    return [=](double timestep) {
        self->set_timestep(timestep);
    };
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<verlet>()
                    .property("integrate", &wrap_integrate<verlet>)
                    .property("finalize", &wrap_finalize<verlet>)
                    .property("set_timestep", &wrap_set_timestep<verlet>)
                    .property("timestep", &verlet::timestep)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("integrate", &runtime::integrate)
                            .def_readonly("finalize", &runtime::finalize)
                    ]
                    .def_readonly("runtime", &verlet::runtime_)

              , def("verlet", &std::make_shared<verlet
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , double
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_integrators_verlet(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    verlet<3, double>::luaopen(L);
    verlet<2, double>::luaopen(L);
#else
    verlet<3, float>::luaopen(L);
    verlet<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class verlet<3, double>;
template class verlet<2, double>;
#else
template class verlet<3, float>;
template class verlet<2, float>;
#endif

} // namespace mdsim
} // namespace host
} // namespace integrators
} // namespace halmd
