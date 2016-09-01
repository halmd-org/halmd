/*
 * Copyright © 2016 Daniel Kirchner
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

#ifndef HALMD_MDSIM_HOST_PARTICLE_ARRAY_HPP
#define HALMD_MDSIM_HOST_PARTICLE_ARRAY_HPP

#include <typeinfo>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/mdsim/host/particle_group.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template<typename T>
class particle_array_typed;

class particle_array : public std::enable_shared_from_this<particle_array>
{
public:
    /**
     * create particle array of given type
     *
     * @param size number of particles
     * @param update_function optional update function
     * @return shared pointer to the new particle array
     */
    template<typename T>
    static inline std::shared_ptr<particle_array_typed<T>> create(unsigned int size,
      std::function<void()> update_function = std::function<void()>());

    /**
     * cast generic particle array to typed data
     *
     * @return typed particle array
     *
     * throws an exception if the type of the generic object does not match
     */
    template<typename T>
    static inline std::shared_ptr<particle_array_typed<T>> cast(std::shared_ptr<particle_array> const&);

    /**
     * empty virtual destructor
     *
     * ensures that the particle array is correctly destroyed even if stored
     * as (shared) pointer to the base class
     */
    virtual ~particle_array() {}

    /**
     * query data type
     *
     * @return RTTI typeid of the stored data
     */
    virtual std::type_info const& type() const = 0;

    /**
    * query cache observer
    *
    * @return a cache observer reflecting the current state of the cache
    */
    virtual cache<> cache_observer() const = 0;

    /**
     * set data from lua table
     *
     * @param object lua table containing the data
     */
    virtual void set_lua(luaponte::object object) = 0;

    /**
     * convert data to lua table
     *
     * @param L lua state (passed in by luaponte)
     * @return lua table containing a copy of the data
     */
    virtual luaponte::object get_lua(lua_State* L) const = 0;
    virtual luaponte::object get_lua(lua_State* L, std::shared_ptr<particle_group> group) const = 0;
    virtual void set_lua(std::shared_ptr<particle_group> group, luaponte::object table) = 0;
};

template<typename T>
class particle_array_typed : public particle_array
{
public:
    /**
     * typed particle array constructor
     *
     * @param args arguments to be passed down to construct the internal raw_array container
     * @param update_function optional update function
     */
    particle_array_typed(unsigned int size, std::function<void()> update_function = std::function<void()>())
      : data_(size), update_function_(update_function)
    {
        if (!update_function_) {
            update_function_ = [](){};
        }
    }

    /**
     * query data type
     *
     * @return RTTI typeid of the stored data
     */
    virtual std::type_info const& type() const
    {
        return typeid(T);
    }

    /**
     * query cache observer
     *
     * @return a cache observer reflecting the current state of the cache
     */
    virtual cache<> cache_observer() const {
        return cache<>(data_);
    }


    /**
     * obtain non-const reference to the stored data
     *
     * @return non-const reference to the cached raw array
     */
    cache<raw_array<T>>& mutable_data()
    {
        return data_;
    }
    /**
     * obtain const reference to the stored data
     *
     * @return const reference to the cached raw array
     */
    cache<raw_array<T>> const& data() const
    {
        update_function_();
        return data_;
    }

    /**
     * set data with iterator
     */
    template <typename iterator_type>
    iterator_type set_data(iterator_type const& first)
    {
        auto output = make_cache_mutable(mutable_data());
        iterator_type input = first;
        for (auto& value : *output) {
            value = *input++;
        }
        return input;
    }
    /**
     * get data with iterator
     */
    template <typename iterator_type>
    iterator_type get_data(iterator_type const& first) const
    {
        auto const& input = read_cache(data());
        return std::copy(input.begin(), input.end(), first);
    }

    /**
     * set data from lua table
     *
     * @param object lua table containing the data
     */
    virtual void set_lua(luaponte::object table)
    {
        auto output = make_cache_mutable(mutable_data());
        std::size_t i = 1;
        for (auto& value : *output) {
            value = luaponte::object_cast<T>(table[i++]);
        }
    }
    /**
     * convert data to lua table
     *
     * @param L lua state (passed in by luaponte)
     * @return lua table containing a copy of the data
     */
    virtual luaponte::object get_lua(lua_State* L) const
    {
        luaponte::object table = luaponte::newtable(L);
        std::size_t i = 1;
        for (auto const& x : read_cache(data())) {
            table[i++] = boost::cref(x);
        }
        return table;
    }

    /**
     * get data from particle group with iterator
     */
    template <typename iterator_type>
    iterator_type get_data(std::shared_ptr<particle_group> particle_group, iterator_type const& first) const
    {
        auto const& group = read_cache(particle_group->ordered());
        auto const& input = read_cache(data());
        auto it = first;
        for(size_t i : group) {
            *it++ = input[i];
        }
        return it;
    }

    /**
     * set data from particle group with iterator
     */
    template <typename iterator_type>
    iterator_type set_data(std::shared_ptr<particle_group> particle_group, iterator_type const& first)
    {
        auto const& group = read_cache(particle_group->ordered());
        auto output = make_cache_mutable(mutable_data());
        auto it = first;
        for(size_t i : group) {
            (*output)[i] = *it++;
        }
        return it;
    }

    virtual void set_lua(std::shared_ptr<particle_group> particle_group, luaponte::object table) {
        auto const& group = read_cache(particle_group->ordered());
        auto output = make_cache_mutable(mutable_data());
        size_t j = 1;
        for(size_t i : group) {
            (*output)[i] = luaponte::object_cast<T>(table[j++]);
        }
    }

    virtual luaponte::object get_lua(lua_State* L, std::shared_ptr<particle_group> particle_group) const {
        auto self = std::static_pointer_cast<particle_array_typed<T> const>(shared_from_this());
        std::function<std::vector<T>()> fn = [self, particle_group]() -> std::vector<T> {
            std::vector<T> data;
            data.reserve(self->data_->size());
            self->get_data(particle_group, std::back_inserter(data));
            return data;
        };
        luaponte::default_converter<std::function<std::vector<T>()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }
private:
    /** cached raw_array */
    cache<raw_array<T>> data_;
    /** optional update function */
    std::function<void()> update_function_;
};

template<typename T>
inline std::shared_ptr<particle_array_typed<T>> particle_array::create(unsigned int size,
  std::function<void()> update_function)
{
    return std::make_shared<particle_array_typed<T>>(size, update_function);
}

template<typename T>
inline std::shared_ptr<particle_array_typed<T>> particle_array::cast(std::shared_ptr<particle_array> const& ptr)
{
    if(ptr->type() != typeid(T)) {
        throw std::invalid_argument("invalid cast of particle array");
    }
    return std::static_pointer_cast<particle_array_typed<T>>(ptr);
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_ARRAY_HPP */
