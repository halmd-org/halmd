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
#include <typeindex>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/gpu/phase_space.hpp>
#include <halmd/observables/gpu/phase_space_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace observables {
namespace gpu {

phase_space_host_cache::phase_space_host_cache(std::shared_ptr<mdsim::gpu::particle_array const> array) : array_(array) {
}

cuda::host::vector<uint8_t>& phase_space_host_cache::acquire(void) {
    if(!(array_->cache_observer() == cache_observer_)) {
        data_ = array_->get_gpu_data();
        cache_observer_ = array_->cache_observer();
    }
    return data_;
}

template<int dimension, typename scalar_type>
class phase_space_sampler_typed : public phase_space_sampler
{
public:
    typedef host::samples::sample<dimension, scalar_type> sample_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::gpu::particle_array_host_wrapper<typename sample_type::data_type> particle_array_type;

    static std::shared_ptr<phase_space_sampler_typed> create(std::shared_ptr<particle_group_type> group
            , std::shared_ptr<mdsim::gpu::particle_array> array
            , std::size_t nthreads
            , std::map<mdsim::gpu::particle_array*, std::shared_ptr<phase_space_host_cache>>& host_cache)
    {
        return std::make_shared<phase_space_sampler_typed>(group, array, nthreads, host_cache);
    }

    phase_space_sampler_typed(std::shared_ptr<particle_group_type> group
            , std::shared_ptr<mdsim::gpu::particle_array> array
            , std::size_t nthreads
            , std::map<mdsim::gpu::particle_array*, std::shared_ptr<phase_space_host_cache>>& host_cache)
            : particle_group_(group), array_(mdsim::gpu::particle_array::cast_host_wrapper<typename sample_type::data_type>(array)),
              nthreads_(nthreads) {
        auto parent = array_->parent();
        auto it = host_cache.find(parent.get());
        if (it == host_cache.end()) {
            host_cache_ = host_cache[parent.get()] = std::make_shared<phase_space_host_cache>(parent);
        } else {
            host_cache_ = it->second;
        }
    }

    virtual std::shared_ptr<host::samples::sample_base> acquire() {
        if (!host_cache_->up_to_date() || group_observer_ != particle_group_->ordered()) {
            auto const& data = host_cache_->acquire();
            auto const& group = particle_group_->ordered_host_cached();

            sample_ = std::make_shared<sample_type>(group.size());
            auto& sample_data = sample_->data();
            auto offset = array_->offset();
            auto stride = array_->stride();

            size_t tag = 0;
            for (size_t i : group) {
                sample_data[tag++] = *reinterpret_cast<typename sample_type::data_type const*>(&data[offset + i * stride]);
            }

            group_observer_ = particle_group_->ordered();

        }
        return sample_;
    }
    virtual void set(std::shared_ptr<host::samples::sample_base const> sample_) {
        if(sample_->type() != typeid(typename sample_type::data_type)) {
            throw std::runtime_error("invalid sample data type");
        }
        auto sample = std::static_pointer_cast<sample_type const>(sample_);
        auto const& sample_data = sample->data();

        auto &data = host_cache_->acquire();
        auto const& group = particle_group_->ordered_host_cached();
        auto offset = array_->offset();
        auto stride = array_->stride();

        typedef typename sample_type::data_type data_type;

        size_t tag = 0;
        for(size_t i : group) {
            *reinterpret_cast<data_type*>(&data[offset + i * stride]) = sample_data[tag++];
        }
#ifdef USE_VERLET_DSFUN
        if(data.capacity() >= data.size() + nthreads_ * stride) {
            offset += nthreads_ * stride;
            for(size_t i : group) {
                *reinterpret_cast<data_type*>(&data[offset + i * stride]) = data_type(0);
            }
        }
#endif
        array_->set_gpu_data(data);
    }

    virtual luaponte::object acquire_lua(lua_State* L, std::shared_ptr<phase_space_sampler> self) {
        std::function<std::shared_ptr<sample_type const>()> fn = [self]() -> std::shared_ptr<sample_type const> {
            return std::static_pointer_cast<sample_type const>(std::static_pointer_cast<phase_space_sampler_typed>(self)->acquire());
        };
        luaponte::default_converter<std::function<std::shared_ptr<sample_type const>()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }
    virtual luaponte::object data_lua(lua_State* L, std::shared_ptr<phase_space_sampler> self) {
        std::function<typename sample_type::array_type const&()> fn = [self]() -> typename sample_type::array_type const& {
            return std::static_pointer_cast<sample_type const>(std::static_pointer_cast<phase_space_sampler_typed>(self)->acquire())->data();
        };
        luaponte::default_converter<std::function<typename sample_type::array_type const&()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }
    virtual void set_lua(luaponte::object sample) {
        set(luaponte::object_cast<std::shared_ptr<sample_type const>>(sample));
    }
protected:
    std::shared_ptr<particle_group_type> particle_group_;
    std::shared_ptr<particle_array_type> array_;
    std::shared_ptr<phase_space_host_cache> host_cache_;
    std::shared_ptr<sample_type> sample_;
    cache<> group_observer_;
    std::size_t nthreads_;
};

static const std::unordered_map<std::type_index,
        std::function<std::shared_ptr<phase_space_sampler>(std::shared_ptr<mdsim::gpu::particle_group>
                , std::shared_ptr<mdsim::gpu::particle_array>
                , std::size_t
                , std::map<mdsim::gpu::particle_array*, std::shared_ptr<phase_space_host_cache>>&)>>
        phase_space_sampler_typed_create_map = {
        { typeid(float), phase_space_sampler_typed<1, float>::create },
        { typeid(fixed_vector<float, 2>), phase_space_sampler_typed<2, float>::create },
        { typeid(fixed_vector<float, 3>), phase_space_sampler_typed<3, float>::create },
        { typeid(fixed_vector<float, 4>), phase_space_sampler_typed<4, float>::create },

        { typeid(double), phase_space_sampler_typed<1, double>::create },
        { typeid(fixed_vector<double, 2>), phase_space_sampler_typed<2, double>::create },
        { typeid(fixed_vector<double, 3>), phase_space_sampler_typed<3, double>::create },
        { typeid(fixed_vector<double, 4>), phase_space_sampler_typed<4, double>::create },

        { typeid(int), phase_space_sampler_typed<1, int>::create },
        { typeid(fixed_vector<int, 2>), phase_space_sampler_typed<2, int>::create },
        { typeid(fixed_vector<int, 3>), phase_space_sampler_typed<3, int>::create },
        { typeid(fixed_vector<int, 4>), phase_space_sampler_typed<4, int>::create },

        { typeid(unsigned int), phase_space_sampler_typed<1, unsigned int>::create },
        { typeid(fixed_vector<unsigned int, 2>), phase_space_sampler_typed<2, unsigned int>::create },
        { typeid(fixed_vector<unsigned int, 3>), phase_space_sampler_typed<3, unsigned int>::create },
        { typeid(fixed_vector<unsigned int, 4>), phase_space_sampler_typed<4, unsigned int>::create },
};


template<int dimension, typename scalar_type>
class phase_space_sampler_position : public phase_space_sampler_typed<dimension, scalar_type>
{
public:
    typedef host::samples::sample<dimension, scalar_type> sample_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle_array_host_wrapper<typename sample_type::data_type> particle_array_type;

    static std::shared_ptr<phase_space_sampler_position> create(std::shared_ptr<particle_group_type> group
            , std::shared_ptr<box_type const> box
            , std::shared_ptr<mdsim::gpu::particle_array> position_array
            , std::shared_ptr<mdsim::gpu::particle_array> image_array
            , cuda::config const& dim
            , std::map<mdsim::gpu::particle_array*, std::shared_ptr<phase_space_host_cache>>& host_cache)
    {
        return std::make_shared<phase_space_sampler_position>(group, box, position_array, image_array, dim, host_cache);
    }

    phase_space_sampler_position(std::shared_ptr<particle_group_type> group
            , std::shared_ptr<box_type const> box
            , std::shared_ptr<mdsim::gpu::particle_array> position_array
            , std::shared_ptr<mdsim::gpu::particle_array> image_array
            , cuda::config const& dim
            , std::map<mdsim::gpu::particle_array*, std::shared_ptr<phase_space_host_cache>>& host_cache)
            : phase_space_sampler_typed<dimension, scalar_type>(group, position_array, dim.threads(), host_cache), box_(box),
              image_array_(mdsim::gpu::particle_array::cast_host_wrapper<typename sample_type::data_type>(image_array)),
              dim_(dim) {
        auto image_parent = image_array_->parent();
        auto it = host_cache.find(image_parent.get());
        if (it == host_cache.end()) {
            image_host_cache_ = host_cache[image_parent.get()] = std::make_shared<phase_space_host_cache>(image_parent);
        } else {
            image_host_cache_ = it->second;
        }
    }

    virtual std::shared_ptr<host::samples::sample_base> acquire() {
        if (!(this->group_observer_ == this->particle_group_->ordered())
            || !this->host_cache_->up_to_date()
            || !image_host_cache_->up_to_date()) {
            auto const& group = this->particle_group_->ordered_host_cached();
            auto const& particle_position = this->host_cache_->acquire();
            auto const& particle_image = image_host_cache_->acquire();
            this->group_observer_ = this->particle_group_->ordered();

            this->sample_ = std::make_shared<sample_type>(group.size());

            auto& sample_position = this->sample_->data();
            auto position_offset = this->array_->offset();
            auto position_stride = this->array_->stride();
            auto image_offset = image_array_->offset();
            auto image_stride = image_array_->stride();

            size_t tag = 0;
            for (size_t i : group) {
                auto& r = sample_position[tag++];
                r = *reinterpret_cast<typename sample_type::data_type const*>(&particle_position[position_offset + i * position_stride]);
                box_->extend_periodic(r, *reinterpret_cast<typename sample_type::data_type const*>(&particle_image[image_offset + i * image_stride]));
            }
        }
        return this->sample_;
    }
    virtual void set(std::shared_ptr<host::samples::sample_base const> sample) {
        // set position data the same as any other kind of data
        phase_space_sampler_typed<dimension, scalar_type>::set(sample);

        // reduce positions on GPU
        typedef typename mdsim::type_traits<dimension, scalar_type>::gpu::coalesced_vector_type gpu_vector_type;
        auto position = make_cache_mutable(mdsim::gpu::particle_array::cast_gpu<float4>(this->array_->parent())->mutable_data());
        auto image = make_cache_mutable(mdsim::gpu::particle_array::cast_gpu<gpu_vector_type>(this->image_array_->parent())->mutable_data());
        auto const& group = read_cache(this->particle_group_->ordered());
        try {
            phase_space_wrapper<dimension>::kernel.r.bind(*position);
            cuda::configure(dim_.grid, dim_.block);
            phase_space_wrapper<dimension>::kernel.reduce_periodic(
                    &*group.begin()
                    , &*position->begin()
                    , &*image->begin()
                    , static_cast<fixed_vector<scalar_type, dimension>>(box_->length())
                    , group.size()
            );
        }
        catch (cuda::error const&)
        {
            LOG_ERROR("failed to reduce particle positions on GPU");
            throw;
        }
    }

private:
    std::shared_ptr<box_type const> box_;
    std::shared_ptr<particle_array_type> image_array_;
    std::shared_ptr<phase_space_host_cache> image_host_cache_;
    cuda::config const dim_;
};

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
phase_space<host_sample<dimension, float_type> >::phase_space(
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
{}

/**
 * Get phase space sampler implementation.
 */
template <int dimension, typename float_type>
std::shared_ptr<phase_space_sampler>
phase_space<host_sample<dimension, float_type> >::get_sampler(std::string const& name)
{
    auto it = samplers_.find(name);

    if (it != samplers_.end()) {
        return it->second;
    } else {
        auto array = particle_->get_array(name);
        if(!name.compare("position")) {
            return (samplers_[name] = phase_space_sampler_position<dimension, float_type>::create
                    (particle_group_, box_, array, particle_->get_array("image"), particle_->dim, host_cache_));
        } else {
            if (array->gpu()) {
                // TODO
                throw std::runtime_error("gpu sampling not yet implemented");
            } else {
                auto it = phase_space_sampler_typed_create_map.find(array->type());
                if(it == phase_space_sampler_typed_create_map.end()) {
                    throw std::runtime_error("invalid sample type");
                }
                return (samplers_[name] = it->second(particle_group_, array, particle_->dim.threads(), host_cache_));
            }
        }
    }
}

template <typename phase_space_type>
static luaponte::object wrap_acquire_host(lua_State* L, std::shared_ptr<phase_space_type> self, std::string const& name)
{
    auto sampler = self->get_sampler(name);
    return sampler->acquire_lua(L, sampler);
}

template <typename phase_space_type>
static void wrap_set(std::shared_ptr<phase_space_type> self, std::string const& name, luaponte::object sample)
{
    self->get_sampler(name)->set_lua(sample);
}

template <typename phase_space_type>
static luaponte::object wrap_data(lua_State* L, std::shared_ptr<phase_space_type> self, std::string const& name)
{
    auto sampler = self->get_sampler(name);
    return sampler->data_lua(L, sampler);
}

template <int dimension, typename float_type>
void phase_space<host_sample<dimension, float_type> >::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space>()
                .def("acquire", &wrap_acquire_host<phase_space>)
                .def("data", &wrap_data<phase_space>)
                .property("dimension", &wrap_dimension<phase_space>)
                .def("set", &wrap_set<phase_space>)
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
    phase_space<host_sample<3, float> >::luaopen(L);
    phase_space<host_sample<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class phase_space<gpu::samples::phase_space<3, float> >;
template class phase_space<gpu::samples::phase_space<2, float> >;
template class phase_space<host_sample<3, float> >;
template class phase_space<host_sample<2, float> >;

} // namespace observables
} // namespace gpu
} // namespace halmd
