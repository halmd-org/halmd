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

namespace halmd {
namespace mdsim {
namespace host {

template<typename T>
class particle_array_typed;

class particle_array
{
public:
    /**
     * create particle array of given type
     *
     * @param args arguments to be forwarded for the construction of the data container
     * @return shared pointer to the new particle array
     */
    template<typename T, typename... Args>
    static inline std::shared_ptr<particle_array_typed<T>> create(Args&&... args);

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
    virtual luaponte::object get_lua(lua_State* L) = 0;

    /**
     * register an on_get handler
     *
     * @param slot function to register for the on_get signal
     * @return signal connection
     *
     * the signal is emitted whenever a non-mutable version of the data is queried
     * and can be used to update out-dated data just in time
     */
    connection on_get(halmd::signal<void()>::slot_function_type const& slot)
    {
        return on_get_.connect(slot);
    }

protected:
    /** on_get signal */
    halmd::signal<void()> on_get_;
};

template<typename T>
class particle_array_typed : public particle_array
{
public:
    /**
     * typed particle array constructor
     *
     * @param args arguments to be passed down to construct the internal raw_array container
     */
    template <typename... Args>
    particle_array_typed(Args&&... args)
    : data_(std::forward<Args>(args)...)
    {
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
        on_get_();
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
    virtual luaponte::object get_lua(lua_State* L)
    {
        on_get_();
        luaponte::object table = luaponte::newtable(L);
        std::size_t i = 1;
        for (auto const& x : read_cache(data())) {
            table[i++] = boost::cref(x);
        }
        return table;
    }
private:
    /** cached raw_array */
    cache<raw_array<T>> data_;
};

template<typename T, typename... Args>
inline std::shared_ptr<particle_array_typed<T>> particle_array::create(Args&& ...args)
{
    return std::make_shared<particle_array_typed<T>>(args...);
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
