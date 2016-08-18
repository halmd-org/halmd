/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <exception>
#include <functional>
#include <memory>

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
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<particle_group_type> particle_group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<clock_type const> clock
  , std::shared_ptr<logger> logger
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
std::shared_ptr<gpu::samples::phase_space<dimension, float_type> const>
phase_space<gpu::samples::phase_space<dimension, float_type> >::acquire()
{
    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    group_array_type const& group = read_cache(particle_group_->ordered());
    position_array_type const& position = read_cache(particle_->position());
    image_array_type const& image = read_cache(particle_->image());
    velocity_array_type const& velocity = read_cache(particle_->velocity());

    scoped_timer_type timer(runtime_.acquire);

    LOG_TRACE("acquire GPU sample");

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = std::make_shared<sample_type>(group.size(), clock_->step());
    }

    phase_space_wrapper<dimension>::kernel.r.bind(position);
    phase_space_wrapper<dimension>::kernel.image.bind(image);
    phase_space_wrapper<dimension>::kernel.v.bind(velocity);

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    phase_space_wrapper<dimension>::kernel.sample(
        &*group.begin()
      , sample_->position()
      , sample_->velocity()
      , static_cast<vector_type>(box_->length())
      , group.size()
    );

    return sample_;
}

template <typename phase_space_type>
static std::function<std::shared_ptr<typename phase_space_type::sample_type const> ()>
wrap_acquire(std::shared_ptr<phase_space_type> self)
{
    return [=]() {
        return self->acquire();
    };
}

template <typename phase_space_type>
static int wrap_dimension(phase_space_type const&)
{
    return phase_space_type::particle_type::vector_type::static_size;
}

template <int dimension, typename float_type>
void phase_space<gpu::samples::phase_space<dimension, float_type> >::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<phase_space>()
                    .property("acquire", &wrap_acquire<phase_space>)
                    .property("dimension", &wrap_dimension<phase_space>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                            .def_readonly("reset", &runtime::reset)
                    ]
                    .def_readonly("runtime", &phase_space::runtime_)

              , def("phase_space", &std::make_shared<phase_space
                   , std::shared_ptr<particle_type>
                   , std::shared_ptr<particle_group_type>
                   , std::shared_ptr<box_type const>
                   , std::shared_ptr<clock_type const>
                   , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

template <int dimension, typename float_type>
phase_space<host::samples::phase_space<dimension, float_type> >::phase_space(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<particle_group_type> particle_group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<clock_type const> clock
  , std::shared_ptr<logger> logger
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
  , threads_(particle_->dim.threads_per_block()) {}

template<int dimension, typename float_type>
typename phase_space<host::samples::phase_space<dimension, float_type>>::group_array_type const&
phase_space<host::samples::phase_space<dimension, float_type> >::read_group_cache_()
{
    // if the indices changed invalidate everything must be renewed
    if (!(group_observer_ == particle_group_->ordered())) {
        position_observer_ = cache<>();
        image_observer_ = cache<>();
        velocity_observer_ = cache<>();
    }

    group_array_type const& group = read_cache(particle_group_->ordered());
    group_observer_ = particle_group_->ordered();

    return group;
}

/**
 * Sample position and species.
 */
template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::acquire_position_species_()
{
    group_array_type const& group = read_group_cache_();
    if (!(position_observer_ == particle_->position()) || !(image_observer_ == particle_->image())) {
        {
            position_array_type const& position = read_cache(particle_->position());
            image_array_type const& image = read_cache(particle_->image());

            cuda::copy(position.begin(), position.end(), h_r_.begin());
            cuda::copy(image.begin(), image.end(), h_image_.begin());
        }

        // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
        // to hold a previous copy of the sample
        position_ = std::make_shared<position_sample_type>(group.size());
        species_ = std::make_shared<species_sample_type>(group.size());

        assert(group.size() == position_->data().size());
        assert(group.size() == species_->data().size());

        typename sample_type::position_array_type& position = position_->data();
        typename sample_type::species_array_type& species = species_->data();

        // copy particle data using reverse tags as on the GPU
        auto const& h_group = particle_group_->ordered_host_cached();

        std::size_t tag = 0;
        for (std::size_t i : h_group) {
            assert(i < h_r_.size());

            // periodically extended particle position
            vector_type r;
            unsigned int type;
            tie(r, type) <<= h_r_[i];
            box_->extend_periodic(r, static_cast<vector_type>(h_image_[i]));

            position[tag] = r;
            species[tag] = type;
            ++tag;
        }

        position_observer_ = particle_->position();
        image_observer_ = particle_->image();
    }
}

/**
 * Sample velocity and mass.
 */
template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::acquire_velocity_mass_()
{
    group_array_type const& group = read_group_cache_();
    if (!(velocity_observer_ == particle_->velocity())) {
        {
            velocity_array_type const& velocity = read_cache(particle_->velocity());
            cuda::copy(velocity.begin(), velocity.end(), h_v_.begin());
        }

        // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
        // to hold a previous copy of the sample
        velocity_ = std::make_shared<velocity_sample_type>(group.size());
        mass_ = std::make_shared<mass_sample_type>(group.size());

        assert(group.size() == velocity_->data().size());
        assert(group.size() == mass_->data().size());

        typename sample_type::velocity_array_type& velocity = velocity_->data();
        typename sample_type::mass_array_type& mass = mass_->data();

        // copy particle data using reverse tags as on the GPU
        auto const& h_group = particle_group_->ordered_host_cached();

        std::size_t tag = 0;
        for (std::size_t i : h_group) {
            assert(i < h_v_.size());
            tie(velocity[tag], mass[tag]) <<= h_v_[i];
            ++tag;
        }

        velocity_observer_ = particle_->velocity();
    }
}

/**
 * Sample position
 */
template <int dimension, typename float_type>
std::shared_ptr<host::samples::sample<dimension, float_type> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire_position()
{
    acquire_position_species_();
    return position_;
}

/**
 * Sample velocity
 */
template <int dimension, typename float_type>
std::shared_ptr<host::samples::sample<dimension, float_type> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire_velocity()
{
    acquire_velocity_mass_();
    return velocity_;
}

/**
 * Sample species
 */
template <int dimension, typename float_type>
std::shared_ptr<host::samples::sample<1, unsigned int> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire_species()
{
    acquire_position_species_();
    return species_;
}

/**
 * Sample mass
 */
template <int dimension, typename float_type>
std::shared_ptr<host::samples::sample<1, float_type> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire_mass()
{
    acquire_velocity_mass_();
    return mass_;
}

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
std::shared_ptr<host::samples::phase_space<dimension, float_type> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire()
{
    // TODO: timing
    acquire_position_species_();
    acquire_velocity_mass_();
    return std::make_shared<host::samples::phase_space<dimension, float_type>>(position_, velocity_, species_, mass_, clock_->step());
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type>>::set_position(typename sample_type::position_array_type const& input) {
    particle_->template set_data<typename particle_type::position_type>("position", particle_group_, input.begin());
    group_array_type const& group = read_cache(particle_group_->ordered());
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());

    try {
        phase_space_wrapper<dimension>::kernel.r.bind(*position);
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        phase_space_wrapper<dimension>::kernel.reduce_periodic(
                &*group.begin()
                , &*position->begin()
                , &*image->begin()
                , static_cast<vector_type>(box_->length())
                , group.size()
        );
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("failed to reduce particle positions on GPU");
        throw;
    }
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type>>::set_species(typename sample_type::species_array_type const& species) {
    particle_->template set_data<typename particle_type::species_type>("species", particle_group_, species.begin());
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type>>::set_mass(typename sample_type::mass_array_type const& mass) {
    particle_->template set_data<typename particle_type::mass_type>("mass", particle_group_, mass.begin());
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type>>::set_velocity(typename sample_type::velocity_array_type const& velocity) {
    particle_->template set_data<typename particle_type::velocity_type>("velocity", particle_group_, velocity.begin());
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::set(std::shared_ptr<sample_type const> sample)
{
    set_position(sample->position());
    set_species(sample->species());
    set_velocity(sample->velocity());
    set_mass(sample->mass());
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::position_array_type const& ()>
wrap_position(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::position_array_type const& {
        return self->acquire()->position();
    };
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::velocity_array_type const& ()>
wrap_velocity(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::velocity_array_type const& {
        return self->acquire()->velocity();
    };
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::species_array_type const& ()>
wrap_species(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::species_array_type const& {
        return self->acquire()->species();
    };
}

template <typename phase_space_type>
static std::function<typename phase_space_type::sample_type::mass_array_type const& ()>
wrap_mass(std::shared_ptr<phase_space_type> self)
{
    return [=]() -> typename phase_space_type::sample_type::mass_array_type const& {
        return self->acquire()->mass();
    };
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space>()
                .property("acquire", &wrap_acquire<phase_space>)
                .property("position", &wrap_position<phase_space>)
                .property("velocity", &wrap_velocity<phase_space>)
                .property("species", &wrap_species<phase_space>)
                .property("mass", &wrap_mass<phase_space>)
                .property("dimension", &wrap_dimension<phase_space>)
                .def("set", &phase_space::set)
                .def("set_position", &phase_space::set_position)
                .def("set_velocity", &phase_space::set_velocity)
                .def("set_mass", &phase_space::set_mass)
                .def("set_species", &phase_space::set_species)
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
                def("phase_space", &std::make_shared<phase_space
                   , std::shared_ptr<particle_type>
                   , std::shared_ptr<particle_group_type>
                   , std::shared_ptr<box_type const>
                   , std::shared_ptr<clock_type const>
                   , std::shared_ptr<logger>
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
