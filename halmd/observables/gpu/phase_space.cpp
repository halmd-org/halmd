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
  , std::shared_ptr<logger_type> logger
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

    cache_proxy<group_array_type const> group = particle_group_->ordered();
    cache_proxy<position_array_type const> position = particle_->position();
    cache_proxy<image_array_type const> image = particle_->image();
    cache_proxy<velocity_array_type const> velocity = particle_->velocity();

    scoped_timer_type timer(runtime_.acquire);

    LOG_TRACE("acquire GPU sample");

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = std::make_shared<sample_type>(group->size(), clock_->step());
    }

    phase_space_wrapper<dimension>::kernel.r.bind(*position);
    phase_space_wrapper<dimension>::kernel.image.bind(*image);
    phase_space_wrapper<dimension>::kernel.v.bind(*velocity);

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    phase_space_wrapper<dimension>::kernel.sample(
        &*group->begin()
      , sample_->position()
      , sample_->velocity()
      , static_cast<vector_type>(box_->length())
      , group->size()
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
                   , std::shared_ptr<logger_type>
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
  , std::shared_ptr<logger_type> logger
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

/**
 * Sample phase_space
 */
template <int dimension, typename float_type>
std::shared_ptr<host::samples::phase_space<dimension, float_type> const>
phase_space<host::samples::phase_space<dimension, float_type> >::acquire()
{
    if (sample_ && sample_->step() == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return sample_;
    }

    cache_proxy<group_array_type const> group = particle_group_->ordered();

    scoped_timer_type timer(runtime_.acquire);

    LOG_TRACE("acquire host sample");

    try {
        cache_proxy<position_array_type const> position = particle_->position();
        cache_proxy<image_array_type const> image = particle_->image();
        cache_proxy<velocity_array_type const> velocity = particle_->velocity();

        cuda::copy(position->begin(), position->end(), h_r_.begin());
        cuda::copy(image->begin(), image->end(), h_image_.begin());
        cuda::copy(velocity->begin(), velocity->end(), h_v_.begin());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy from GPU to host");
        throw;
    }

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    {
        scoped_timer_type timer(runtime_.reset);
        sample_ = std::make_shared<sample_type>(group->size(), clock_->step());
    }

    assert(group->size() == sample_->position().size());
    assert(group->size() == sample_->velocity().size());
    assert(group->size() == sample_->species().size());
    assert(group->size() == sample_->mass().size());

    typename sample_type::position_array_type& position = sample_->position();
    typename sample_type::velocity_array_type& velocity = sample_->velocity();
    typename sample_type::species_array_type& species = sample_->species();
    typename sample_type::mass_array_type& mass = sample_->mass();

    // copy particle data using reverse tags as on the GPU
    cuda::host::vector<unsigned int> h_group(group->size());
    cuda::copy(group->begin(), group->end(), h_group.begin());

    std::size_t tag = 0;
    for (std::size_t i : h_group) {
        assert(i < h_r_.size());

        // periodically extended particle position
        vector_type r;
        unsigned int type;
        tie(r, type) <<= h_r_[i];
        box_->extend_periodic(r, static_cast<vector_type>(h_image_[i]));

        position[tag] = r;
        tie(velocity[tag], mass[tag]) <<= h_v_[i];
        species[tag] = type;
        ++tag;
    }

    return sample_;
}

template <int dimension, typename float_type>
void phase_space<host::samples::phase_space<dimension, float_type> >::set(std::shared_ptr<sample_type const> sample)
{
    cache_proxy<group_array_type const> group = particle_group_->ordered();
    cache_proxy<position_array_type> position = particle_->position();
    cache_proxy<image_array_type> image = particle_->image();
    cache_proxy<velocity_array_type> velocity = particle_->velocity();

    scoped_timer_type timer(runtime_.set);

    // allocate additional memory for double-single precision
    h_r_.reserve(position->capacity());
    h_v_.reserve(velocity->capacity());

    // copy particle arrays from GPU to host
    cuda::copy(position->begin(), position->begin() + position->capacity(), h_r_.begin());
    cuda::copy(velocity->begin(), velocity->begin() + velocity->capacity(), h_v_.begin());

    // assign particle coordinates and types
    typename particle_type::species_type const nspecies = particle_->nspecies();
    typename sample_type::position_array_type const& sample_position = sample->position();
    typename sample_type::velocity_array_type const& sample_velocity = sample->velocity();
    typename sample_type::species_array_type const& sample_species = sample->species();
    typename sample_type::mass_array_type const& sample_mass = sample->mass();

    assert(sample_position.size() >= group->size());
    assert(sample_velocity.size() >= group->size());
    assert(sample_species.size() >= group->size());

    // copy particle data using reverse tags as on the GPU
    cuda::host::vector<unsigned int> h_group(group->size());
    cuda::copy(group->begin(), group->end(), h_group.begin());

    std::size_t tag = 0;
    std::size_t nthreads = particle_->dim.threads();
    for (std::size_t i : h_group) {
        assert(i < h_r_.size());
        unsigned int species = sample_species[tag];
        if (species >= nspecies) {
            throw std::invalid_argument("invalid species");
        }
        h_r_[i] <<= tie(sample_position[tag], species);
        h_v_[i] <<= tie(sample_velocity[tag], sample_mass[tag]);
#ifdef USE_VERLET_DSFUN
        h_r_[i + nthreads] = vector_type(0);
        h_v_[i + nthreads] = vector_type(0);
#endif
        ++tag;
    }

    try {
        cuda::copy(h_r_.begin(), h_r_.begin() + h_r_.capacity(), position->begin());
        cuda::copy(h_v_.begin(), h_v_.begin() + h_v_.capacity(), velocity->begin());
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("failed to copy particles to GPU");
        throw;
    }

    // shift particle positions to range (-L/2, L/2)
    try {
        phase_space_wrapper<dimension>::kernel.r.bind(*position);
        cuda::configure(particle_->dim.grid, particle_->dim.block);
        phase_space_wrapper<dimension>::kernel.reduce_periodic(
            &*group->begin()
          , &*position->begin()
          , &*image->begin()
          , static_cast<vector_type>(box_->length())
          , group->size()
        );
    }
    catch (cuda::error const&)
    {
        LOG_ERROR("failed to reduce particle positions on GPU");
        throw;
    }
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

/**
 * Translate species from internal (0-based) to external (1-based) representation.
 *
 * This function returns a const reference to a species array that is
 * stored in the functor, and passed by reference to this function.
 */
template <typename phase_space_type>
static std::function<std::vector<typename phase_space_type::sample_type::species_array_type::value_type> const& ()>
wrap_species(std::shared_ptr<phase_space_type> self)
{
    typedef typename phase_space_type::sample_type sample_type;
    typedef typename sample_type::species_array_type::value_type species_type;

    std::shared_ptr<std::vector<species_type>> species = std::make_shared<std::vector<species_type>>();
    return [=]() -> std::vector<species_type> const& {
        std::shared_ptr<sample_type const> sample = self->acquire();
        species->resize(sample->species().size());
        std::transform(
            sample->species().begin()
          , sample->species().end()
          , species->begin()
          , [](species_type s) {
                return s + 1;
            }
        );
        return *species;
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

template <typename phase_space_type>
static std::function<void (std::shared_ptr<typename phase_space_type::sample_type const>)>
wrap_set(std::shared_ptr<phase_space_type> self)
{
    return [=](std::shared_ptr<typename phase_space_type::sample_type const> sample) {
        self->set(sample);
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
                .property("set", &wrap_set<phase_space>)
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
                   , std::shared_ptr<logger_type>
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
