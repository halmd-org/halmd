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
#include <functional>
#include <memory>
#include <typeindex>

#include <halmd/observables/host/phase_space.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace observables {
namespace host {

/**
 * phase space sampler implementation for typed samples
 *
 * copies the data to a host sample and provides the actual implementation
 * for the sample related interface of phase_space
 */
template<int dimension, typename scalar_type>
class phase_space_sampler_typed : public phase_space_sampler
{
public:
    typedef samples::sample<dimension, scalar_type> sample_type;
    typedef mdsim::host::particle_group particle_group_type;
    typedef mdsim::host::particle_array_typed<typename sample_type::data_type> particle_array_type;

    static std::shared_ptr<phase_space_sampler_typed> create(std::shared_ptr<particle_group_type> group
                                               , std::shared_ptr<mdsim::host::particle_array> array)
    {
        return std::make_shared<phase_space_sampler_typed>(group, array);
    }

    phase_space_sampler_typed(std::shared_ptr<particle_group_type> group
            , std::shared_ptr<mdsim::host::particle_array> array)
            : particle_group_(group), array_(mdsim::host::particle_array::cast<typename sample_type::data_type>(array)) {
    }

    /**
     * acquire a sample
     *
     * copies the data to a new sample if the data is not up to date,
     * returns the stored sample otherwise
     */
    virtual std::shared_ptr<sample_base> acquire() {
        if (!(group_observer_ == particle_group_->ordered()) || !(array_observer_ == array_->cache_observer())) {
            auto const& group = read_cache(particle_group_->ordered());
            group_observer_ = particle_group_->ordered();

            auto const& data = read_cache(array_->data());

            sample_ = std::make_shared<sample_type>(group.size());

            auto& sample_data = sample_->data();
            // copy velocities using index map
            std::size_t tag = 0;
            for (std::size_t i : group) {
                sample_data[tag] = data[i];
                ++tag;
            }
            array_observer_ = array_->cache_observer();
        }
        return sample_;
    }
    /**
     * copies the data of a sample to the particle array
     */
    virtual void set(std::shared_ptr<sample_base const> sample) {
        if(sample->type() != typeid(typename sample_type::data_type)) {
            throw std::runtime_error("invalid sample data type");
        }
        auto const& group = read_cache(particle_group_->ordered());
        auto data = make_cache_mutable(array_->mutable_data());
        auto const& sample_data = std::static_pointer_cast<sample_type const>(sample)->data();

        std::size_t tag = 0;
        for (std::size_t i : group) {
            (*data)[i] = sample_data[tag];
            ++tag;
        }
    }

    /**
     * returns a lua slot function to be used to acquire a host sample
     */
    virtual luaponte::object acquire_lua(lua_State* L, std::shared_ptr<phase_space_sampler> self) {
        std::function<std::shared_ptr<sample_type const>()> fn = [self]() -> std::shared_ptr<sample_type const> {
            return std::static_pointer_cast<sample_type const>(std::static_pointer_cast<phase_space_sampler_typed>(self)->acquire());
        };
        luaponte::default_converter<std::function<std::shared_ptr<sample_type const>()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }

    /**
     * returns a lua slot function to be used to directly acquire the data of a host sample
     */
    virtual luaponte::object data_lua(lua_State* L, std::shared_ptr<phase_space_sampler> self) {
        std::function<typename sample_type::array_type const&()> fn = [self]() -> typename sample_type::array_type const& {
            return std::static_pointer_cast<sample_type const>(std::static_pointer_cast<phase_space_sampler_typed>(self)->acquire())->data();
        };
        luaponte::default_converter<std::function<typename sample_type::array_type const&()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }

    /**
     * wrapper to export the set member to lua
     */
    virtual void set_lua(luaponte::object sample) {
        set(luaponte::object_cast<std::shared_ptr<sample_type const>>(sample));
    }
protected:
    std::shared_ptr<particle_group_type> particle_group_;
    std::shared_ptr<particle_array_type> array_;
    std::shared_ptr<sample_type> sample_;
    cache<> array_observer_;
    cache<> group_observer_;
};

/**
 * associative map from typeid's to typed phase space sampler create functions
 *
 * used to create the correct phase space sampler based on the typeid of a particle array
 */
static const std::unordered_map<std::type_index,
    std::function<std::shared_ptr<phase_space_sampler>(std::shared_ptr<mdsim::host::particle_group>
                                                     , std::shared_ptr<mdsim::host::particle_array>)>>
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
    typedef samples::sample<dimension, scalar_type> sample_type;
    typedef mdsim::host::particle_group particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle_array_typed<typename sample_type::data_type> particle_array_type;

    static std::shared_ptr<phase_space_sampler_position> create(std::shared_ptr<particle_group_type> group
            , std::shared_ptr<box_type const> box
            , std::shared_ptr<mdsim::host::particle_array> position_array
            , std::shared_ptr<mdsim::host::particle_array> image_array)
    {
        return std::make_shared<phase_space_sampler_position>(group, box, position_array, image_array);
    }

    phase_space_sampler_position(std::shared_ptr<particle_group_type> group
    , std::shared_ptr<box_type const> box
    , std::shared_ptr<mdsim::host::particle_array> position_array
    , std::shared_ptr<mdsim::host::particle_array> image_array)
    : phase_space_sampler_typed<dimension, scalar_type>(group, position_array), box_(box),
      image_array_(mdsim::host::particle_array::cast<typename sample_type::data_type>(image_array)) {
    }

    virtual std::shared_ptr<sample_base> acquire() {
        if (!(this->group_observer_ == this->particle_group_->ordered())
            || !(this->array_observer_ == this->array_->cache_observer())
            || !(image_array_observer_ == image_array_->cache_observer())) {
            auto const& group = read_cache(this->particle_group_->ordered());
            this->group_observer_ = this->particle_group_->ordered();

            auto const& particle_position = read_cache(this->array_->data());
            auto const& particle_image = read_cache(image_array_->data());

            this->sample_ = std::make_shared<sample_type>(group.size());

            auto& sample_position = this->sample_->data();

            // copy and periodically extend positions using index map
            std::size_t tag = 0;
            for (std::size_t i : group) {
                auto& r = sample_position[tag];
                r = particle_position[i];
                box_->extend_periodic(r, particle_image[i]);
                ++tag;
            }

            this->array_observer_ = this->array_->cache_observer();
            image_array_observer_ = image_array_->cache_observer();
        }
        return this->sample_;
    }
    virtual void set(std::shared_ptr<sample_base const> sample) {
        if(sample->type() != typeid(typename sample_type::data_type)) {
            throw std::runtime_error("invalid sample data type");
        }

        auto const& group = read_cache(this->particle_group_->ordered());
        auto particle_position = make_cache_mutable(this->array_->mutable_data());
        auto particle_image = make_cache_mutable(image_array_->mutable_data());

        auto const& sample_position = std::static_pointer_cast<sample_type const>(sample)->data();

        std::size_t tag = 0;
        for (std::size_t i : group) {
            auto& r = (*particle_position)[i];
            r = sample_position[tag];
            auto& image = (*particle_image)[i];
            image = 0;

            // The host implementation of reduce_periodic wraps the position at
            // most once around the box. This is more efficient during the
            // simulation. For setting arbitrary particle positions, however,
            // we must ensure that the final position is actually inside the
            // periodic box.
            typedef typename sample_type::data_type vector_type;
            vector_type shift;
            do {
                image += (shift = box_->reduce_periodic(r));
            } while (shift != vector_type(0));
            ++tag;
        }
    }
private:
    std::shared_ptr<box_type const> box_;
    std::shared_ptr<particle_array_type> image_array_;
    cache<> image_array_observer_;
};

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<particle_group_type> particle_group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , particle_group_(particle_group)
  , box_(box)
  , logger_(logger) {}

/**
 * Get phase space sampler implementation.
 */
template <int dimension, typename float_type>
std::shared_ptr<phase_space_sampler>
phase_space<dimension, float_type>::get_sampler(std::string const& name)
{
    auto it = samplers_.find(name);
    if (it != samplers_.end()) {
        return it->second;
    } else {
        auto array = particle_->get_array(name);
        if(!name.compare("position")) {
            auto image = particle_->get_array("image");
            auto sampler = phase_space_sampler_position<dimension, float_type>::create(particle_group_, box_, array, image);
            return (samplers_[name] = sampler);
        } else {
            auto it = phase_space_sampler_typed_create_map.find(array->type());
            if(it == phase_space_sampler_typed_create_map.end()) {
                throw std::runtime_error("invalid sample type");
            }
            return (samplers_[name] = it->second(particle_group_, array));
        }
    }
}

template <typename phase_space_type>
static luaponte::object wrap_acquire(lua_State* L, std::shared_ptr<phase_space_type> self, std::string const& name)
{
    auto sampler = self->get_sampler(name);
    return sampler->acquire_lua(L, sampler);
}

template <typename phase_space_type>
static luaponte::object wrap_data(lua_State* L, std::shared_ptr<phase_space_type> self, std::string const& name)
{
    auto sampler = self->get_sampler(name);
    return sampler->data_lua(L, sampler);
}

template <typename phase_space_type>
static void wrap_set(std::shared_ptr<phase_space_type> self, std::string const& name, luaponte::object sample)
{
    self->get_sampler(name)->set_lua(sample);
}

template <typename phase_space_type>
static int wrap_dimension(phase_space_type const&)
{
    return phase_space_type::particle_type::vector_type::static_size;
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space>()
                .def("acquire", &wrap_acquire<phase_space>)
                .def("data", &wrap_data<phase_space>)
                .def("set", &wrap_set<phase_space>)
                .property("dimension", &wrap_dimension<phase_space>)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                        .def_readonly("reset", &runtime::reset)
                        .def_readonly("set", &runtime::set)
                ]
                .def_readonly("runtime", &phase_space::runtime_)

          , def("phase_space", &std::make_shared<phase_space
               , std::shared_ptr<particle_type>
               , std::shared_ptr<particle_group_type>
               , std::shared_ptr<box_type const>
               , std::shared_ptr<logger>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    phase_space<3, double>::luaopen(L);
    phase_space<2, double>::luaopen(L);
#else
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#else
template class phase_space<3, float>;
template class phase_space<2, float>;
#endif

} // namespace observables
} // namespace host
} // namespace halmd
