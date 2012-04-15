/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/observables/gpu/phase_space.hpp>
#include <halmd/observables/gpu/phase_space_kernel.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
phase_space<gpu::samples::phase_space<dimension, float_type> >::phase_space(
    shared_ptr<particle_group_type const> particle_group
  , shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_group_(particle_group)
  , particle_(particle_group->particle())
  , box_(box)
  , clock_(clock)
  , logger_(logger)
{
    try {
        cuda::copy(static_cast<vector_type>(box_->length()), phase_space_wrapper<dimension>::kernel.box_length);
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("failed to copy box length to GPU");
        throw;
    }
}

template <int dimension, typename float_type>
phase_space<host::samples::phase_space<dimension, float_type> >::phase_space(
    shared_ptr<particle_group_type> particle_group
  , shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_group_(particle_group)
  , particle_(particle_group->particle())
  , box_(box)
  , clock_(clock)
  , logger_(logger)
  // allocate page-locked host memory
  , h_r_(particle_->nbox)
  , h_image_(particle_->nbox)
  , h_v_(particle_->nbox)
{
}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
shared_ptr<gpu::samples::phase_space<dimension, float_type> const>
phase_space<gpu::samples::phase_space<dimension, float_type> >::acquire()
{
    scoped_timer_type timer(runtime_.acquire);

    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    LOG_TRACE("acquire GPU sample");

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = make_shared<sample_type>(particle_group_->size(), clock_->step());
    }

    phase_space_wrapper<dimension>::kernel.r.bind(particle_->g_r);
    phase_space_wrapper<dimension>::kernel.image.bind(particle_->g_image);
    phase_space_wrapper<dimension>::kernel.v.bind(particle_->g_v);

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    phase_space_wrapper<dimension>::kernel.sample(
        particle_group_->g_map()
      , sample_->position()
      , sample_->velocity()
      , particle_group_->size()
    );

    return sample_;
}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
shared_ptr<host::samples::phase_space<dimension, float_type> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire()
{
    scoped_timer_type timer(runtime_.acquire);

    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    LOG_TRACE("acquire host sample");

    try {
        cuda::copy(particle_->g_r, h_r_);
        cuda::copy(particle_->g_image, h_image_);
        cuda::copy(particle_->g_v, h_v_);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy from GPU to host");
        throw;
    }

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = make_shared<sample_type>(particle_group_->size(), clock_->step());
    }

    assert(particle_group_->size() == sample_->position().size());
    assert(particle_group_->size() == sample_->velocity().size());
    assert(particle_group_->size() == sample_->species().size());

    typename sample_type::position_array_type& position = sample_->position();
    typename sample_type::velocity_array_type& velocity = sample_->velocity();
    typename sample_type::species_array_type& species = sample_->species();

    // copy particle data using reverse tags as on the GPU
    unsigned int const* map = particle_group_->h_map();
    for (unsigned int tag = 0; tag < particle_group_->size(); ++tag) {

        unsigned int idx = map[tag];
        assert(idx < h_r_.size());

        using mdsim::gpu::particle_kernel::untagged;

        // periodically extended particle position
        vector_type r;
        unsigned int type;
        tie(r, type) = untagged<vector_type>(h_r_[idx]);
        box_->extend_periodic(r, static_cast<vector_type>(h_image_[idx]));

        position[tag] = r;
        velocity[tag] = static_cast<vector_type>(h_v_[idx]);
        species[tag] = type;
    }

    return sample_;
}

template <typename phase_space_type, typename sample_type>
static function<shared_ptr<sample_type const> ()>
wrap_acquire(shared_ptr<phase_space_type> phase_space)
{
    return bind(&phase_space_type::acquire, phase_space);
}

template <typename phase_space_type>
static int wrap_dimension(phase_space_type const&)
{
    return phase_space_type::particle_type::vector_type::static_size;
}

template <int dimension, typename float_type>
void phase_space<gpu::samples::phase_space<dimension, float_type> >::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name(demangled_name<phase_space>());
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space>(class_name.c_str())
                .property("acquire", &wrap_acquire<phase_space, sample_type>)
                .property("dimension", &wrap_dimension<phase_space>)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                        .def_readonly("reset", &runtime::reset)
                ]
                .def_readonly("runtime", &phase_space::runtime_)

          , namespace_("gpu")
            [
                def("phase_space", &make_shared<phase_space
                   , shared_ptr<particle_group_type const>
                   , shared_ptr<box_type const>
                   , shared_ptr<clock_type const>
                   , shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name(demangled_name<phase_space>());
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space>(class_name.c_str())
                .property("acquire", &wrap_acquire<phase_space, sample_type>)
                .property("dimension", &wrap_dimension<phase_space>)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                        .def_readonly("reset", &runtime::reset)
                ]
                .def_readonly("runtime", &phase_space::runtime_)

          , namespace_("host")
            [
                def("phase_space", &make_shared<phase_space
                   , shared_ptr<particle_group_type>
                   , shared_ptr<box_type const>
                   , shared_ptr<clock_type const>
                   , shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_phase_space(lua_State* L)
{
    phase_space<gpu::samples::phase_space<3, float> >::luaopen(L);
    phase_space<gpu::samples::phase_space<2, float> >::luaopen(L);
    phase_space<host::samples::phase_space<3, float> >::luaopen(L);
    phase_space<host::samples::phase_space<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class phase_space<gpu::samples::phase_space<3, float> >;
template class phase_space<gpu::samples::phase_space<2, float> >;
template class phase_space<host::samples::phase_space<3, float> >;
template class phase_space<host::samples::phase_space<2, float> >;

} // namespace observables
} // namespace gpu
} // namespace halmd
