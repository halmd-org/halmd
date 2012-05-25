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

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/observables/gpu/phase_space.hpp>
#include <halmd/observables/gpu/phase_space_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
phase_space<gpu::samples::phase_space<dimension, float_type> >::phase_space(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<particle_group_type const> particle_group
  , boost::shared_ptr<box_type const> box
  , boost::shared_ptr<clock_type const> clock
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , particle_group_(particle_group)
  , box_(box)
  , clock_(clock)
  , logger_(logger) {}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
boost::shared_ptr<gpu::samples::phase_space<dimension, float_type> const>
phase_space<gpu::samples::phase_space<dimension, float_type> >::acquire()
{
    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    scoped_timer_type timer(runtime_.acquire);

    LOG_TRACE("acquire GPU sample");

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = boost::make_shared<sample_type>(particle_group_->size(), clock_->step());
    }

    phase_space_wrapper<dimension>::kernel.r.bind(particle_->position());
    phase_space_wrapper<dimension>::kernel.image.bind(particle_->image());
    phase_space_wrapper<dimension>::kernel.v.bind(particle_->velocity());

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    phase_space_wrapper<dimension>::kernel.sample(
        particle_group_->begin()
      , sample_->position()
      , sample_->velocity()
      , static_cast<vector_type>(box_->length())
      , particle_group_->size()
    );

    return sample_;
}

template <typename phase_space_type, typename sample_type>
static boost::function<boost::shared_ptr<sample_type const> ()>
wrap_acquire(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&phase_space_type::acquire, phase_space);
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
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<phase_space>()
                    .property("acquire", &wrap_acquire<phase_space, sample_type>)
                    .property("dimension", &wrap_dimension<phase_space>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                            .def_readonly("reset", &runtime::reset)
                    ]
                    .def_readonly("runtime", &phase_space::runtime_)

              , def("phase_space", &boost::make_shared<phase_space
                   , boost::shared_ptr<particle_type>
                   , boost::shared_ptr<particle_group_type const>
                   , boost::shared_ptr<box_type const>
                   , boost::shared_ptr<clock_type const>
                   , boost::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

template <int dimension, typename float_type>
phase_space<host::samples::phase_space<dimension, float_type> >::phase_space(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<particle_group_type const> particle_group
  , boost::shared_ptr<box_type const> box
  , boost::shared_ptr<clock_type const> clock
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , particle_group_(particle_group)
  , box_(box)
  , clock_(clock)
  , logger_(logger)
  // allocate page-locked host memory
  , h_r_(particle_->nparticle())
  , h_image_(particle_->nparticle())
  , h_v_(particle_->nparticle())
  , h_group_(particle_group_->size())
  , threads_(particle_->dim.threads_per_block()) {}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
boost::shared_ptr<host::samples::phase_space<dimension, float_type> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire()
{
    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    scoped_timer_type timer(runtime_.acquire);

    LOG_TRACE("acquire host sample");

    try {
        cuda::copy(particle_->position(), h_r_);
        cuda::copy(particle_->image(), h_image_);
        cuda::copy(particle_->velocity(), h_v_);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy from GPU to host");
        throw;
    }

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = boost::make_shared<sample_type>(particle_group_->size(), clock_->step());
    }

    assert(particle_group_->size() == sample_->position().size());
    assert(particle_group_->size() == sample_->velocity().size());
    assert(particle_group_->size() == sample_->species().size());

    typename sample_type::position_array_type& position = sample_->position();
    typename sample_type::velocity_array_type& velocity = sample_->velocity();
    typename sample_type::species_array_type& species = sample_->species();

    // copy particle data using reverse tags as on the GPU
    cuda::configure((particle_group_->size() + threads_ - 1) / threads_, threads_);
    phase_space_wrapper<dimension>::kernel.copy_particle_group(
        particle_group_->begin()
      , h_group_
      , particle_group_->size()
    );
    cuda::thread::synchronize();

    std::size_t tag = 0;
    for (std::size_t i : h_group_) {
        assert(i < h_r_.size());

        // periodically extended particle position
        vector_type r;
        unsigned int type;
        tie(r, type) <<= h_r_[i];
        box_->extend_periodic(r, static_cast<vector_type>(h_image_[i]));

        position[tag] = r;
        velocity[tag] = static_cast<vector_type>(h_v_[i]);
        species[tag] = type;
        ++tag;
    }

    return sample_;
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::set(boost::shared_ptr<sample_type const> sample)
{
    scoped_timer_type timer(runtime_.set);

    // assign particle coordinates and types
    typename sample_type::position_array_type const& sample_position = sample->position();
    typename sample_type::velocity_array_type const& sample_velocity = sample->velocity();
    typename sample_type::species_array_type const& sample_species = sample->species();

    assert(sample_position.size() >= particle_group_->size());
    assert(sample_velocity.size() >= particle_group_->size());
    assert(sample_species.size() >= particle_group_->size());

    // copy particle data using reverse tags as on the GPU
    cuda::configure((particle_group_->size() + 128 - 1) / 128, 128);
    phase_space_wrapper<dimension>::kernel.copy_particle_group(
        particle_group_->begin()
      , h_group_
      , particle_group_->size()
    );
    cuda::thread::synchronize();

    std::size_t tag = 0;
    for (std::size_t i : h_group_) {
        assert(i < h_r_.size());
        h_r_[i] <<= tie(sample_position[tag], sample_species[tag]);
        h_v_[i] = sample_velocity[tag];
        ++tag;
    }

    try {
#ifdef USE_VERLET_DSFUN
        cuda::memset(particle_->position(), 0, particle_->position().capacity());
        cuda::memset(particle_->velocity(), 0, particle_->velocity().capacity());
#endif
        cuda::copy(h_r_, particle_->position());
        cuda::copy(h_v_, particle_->velocity());
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("failed to copy particles to GPU");
        throw;
    }

    // shift particle positions to range (-L/2, L/2)
    try {
        phase_space_wrapper<dimension>::kernel.r.bind(particle_->position());
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        phase_space_wrapper<dimension>::kernel.reduce_periodic(
            particle_group_->begin()
          , particle_->position()
          , particle_->image()
          , static_cast<vector_type>(box_->length())
          , particle_group_->size()
        );
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("failed to reduce particle positions on GPU");
        throw;
    }
}

template <typename phase_space_type>
static typename phase_space_type::sample_type::position_array_type const&
position(boost::shared_ptr<phase_space_type> const& phase_space)
{
    return phase_space->acquire()->position();
}

template <typename phase_space_type>
static boost::function<typename phase_space_type::sample_type::position_array_type const& ()>
wrap_position(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&position<phase_space_type>, phase_space);
}

template <typename phase_space_type>
static typename phase_space_type::sample_type::velocity_array_type const&
velocity(boost::shared_ptr<phase_space_type> const& phase_space)
{
    return phase_space->acquire()->velocity();
}

template <typename phase_space_type>
static boost::function<typename phase_space_type::sample_type::velocity_array_type const& ()>
wrap_velocity(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&velocity<phase_space_type>, phase_space);
}

/**
 * Translate species from internal (0-based) to external (1-based) representation.
 *
 * This function returns a const reference to a species array that is
 * stored in the functor, and passed by reference to this function.
 */
template <typename phase_space_type>
static typename phase_space_type::sample_type::species_array_type const&
species(boost::shared_ptr<phase_space_type> const& phase_space, boost::shared_ptr<typename phase_space_type::sample_type::species_array_type>& species)
{
    boost::shared_ptr<typename phase_space_type::sample_type const> sample = phase_space->acquire();
    species.reset(new typename phase_space_type::sample_type::species_array_type(sample->species().size()));
    std::transform(
        sample->species().begin()
      , sample->species().end()
      , species->begin()
      , [](typename phase_space_type::sample_type::species_array_type::value_type s) {
            return s + 1;
        }
    );
    return *species;
}

template <typename phase_space_type>
static boost::function<typename phase_space_type::sample_type::species_array_type const& ()>
wrap_species(boost::shared_ptr<phase_space_type> phase_space)
{
    return boost::bind(&species<phase_space_type>, phase_space, boost::shared_ptr<typename phase_space_type::sample_type::species_array_type>());
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space>()
                .property("acquire", &wrap_acquire<phase_space, sample_type>)
                .property("position", &wrap_position<phase_space>)
                .property("velocity", &wrap_velocity<phase_space>)
                .property("species", &wrap_species<phase_space>)
                .property("dimension", &wrap_dimension<phase_space>)
                .def("set", &phase_space::set)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                        .def_readonly("reset", &runtime::reset)
                        .def_readonly("set", &runtime::set)
                ]
                .def_readonly("runtime", &phase_space::runtime_)

          , namespace_("host")
            [
                def("phase_space", &boost::make_shared<phase_space
                   , boost::shared_ptr<particle_type>
                   , boost::shared_ptr<particle_group_type const>
                   , boost::shared_ptr<box_type const>
                   , boost::shared_ptr<clock_type const>
                   , boost::shared_ptr<logger_type>
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
