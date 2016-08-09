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

#ifndef HALMD_MDSIM_GPU_PARTICLE_DATA_HPP
#define HALMD_MDSIM_GPU_PARTICLE_DATA_HPP

#include <typeinfo>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

// forward declarations
template<typename T>
class particle_data_typed;
template<typename T>
class particle_data_gpu;
template<typename host_type, typename gpu_type, bool is_tuple = false>
class particle_data_host_wrapper;

/** particle data base class */
class particle_data
{
public:
    /**
     * create gpu particle data of given type
     *
     * @param size number of particles
     * @param update_function optional update function
     * @return shared pointer to the new particle data
     */
    template<typename T>
    static inline std::shared_ptr<particle_data_gpu<T>> create(unsigned int size, std::function<void()> update_function = std::function<void()>());

    /**
     * cast generic particle data to typed data
     *
     * @param ptr generic particle data
     * @return typed particle data
     *
     * throws an exception if the type of the generic object does not match
     */
    template<typename T>
    static inline std::shared_ptr<particle_data_typed<T>> cast(std::shared_ptr<particle_data> const& ptr);

    /**
     * cast generic particle data to typed gpu data
     *
     * @param ptr generic particle data
     * @return typed gpu particle data
     *
     * throws an exception if the type of the generic object does not match
     */
    template<typename T>
    static inline std::shared_ptr<particle_data_gpu<T>> cast_gpu(std::shared_ptr<particle_data> const& ptr);

    /**
     * create host wrapper for particle data of given gpu and host type
     *
     * @param parent gpu particle data
     * @return shared pointer to the particle data wrapper
     *
     * this version returns an actual wrapper object for the case that gpu_type and host_type
     * are actually different
     */
    template<typename host_type, typename gpu_type>
    static inline typename std::enable_if<!std::is_same<host_type, gpu_type>::value,
            std::shared_ptr<particle_data_host_wrapper<host_type, gpu_type>>>::type
    create_host_wrapper(std::shared_ptr<particle_data_gpu<gpu_type>> const& parent);

    /**
     * create host wrapper for particle data of given gpu and host type
     *
     * @param parent gpu particle data
     * @return shared pointer to the particle data wrapper
     *
     * this version directly passes down the parent object for the case that gpu_type and host_type
     * are actually the same (e.g. in the case of scalar values)
     */
    template<typename host_type, typename gpu_type>
    static inline typename std::enable_if<std::is_same<host_type, gpu_type>::value,
            std::shared_ptr<particle_data_gpu<gpu_type>>>::type
    create_host_wrapper(std::shared_ptr<particle_data_gpu<gpu_type>> const& parent);

    /**
     * create host wrapper for a tuple element of packed gpu particle data
    *
    * @param parent gpu particle data
    * @return shared pointer to the particle data wrapper
    */
    template<typename tuple_type, int field, typename gpu_type>
    static inline std::shared_ptr<particle_data_host_wrapper<typename std::tuple_element<field, tuple_type>::type, gpu_type, true>>
    create_packed_wrapper(std::shared_ptr<particle_data_gpu<gpu_type>> const& parent);

    /**
     * empty virtual destructor
     *
     * ensures that the particle data is correctly destroyed even if stored
     * as (shared) pointer to the base class
     */
    virtual ~particle_data(){}

    /**
     * query data type
     *
     * @return RTTI typeid of the stored data
     */
    virtual std::type_info const& type() const = 0;
    /**
     * query gpu flag
     *
     * @return true if the particle data is raw gpu data, false otherwise
     */
    virtual bool gpu() const = 0;

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
};

/** particle data base class with type information */
template<typename T>
class particle_data_typed : public particle_data
{
public:
    /**
     * typed particle data constructor
     *
     * @param stride stride of the underlying gpu data
     * @param offset offset of the typed data within the underlying gpu data
     */
    particle_data_typed(size_t stride, size_t offset) : stride_(stride), offset_(offset) {
    }
    /**
     * query data type
     *
     * @return RTTI typeid of the stored data
     */
    virtual std::type_info const& type() const {
        return typeid(T);
    }

    /**
     * set data with iterator
     */
    template <typename iterator_type>
    iterator_type set_data(iterator_type const& first)
    {
        auto mem = get_memory();
        auto it = first;
        for(size_t i = offset_; i < mem.size(); i += stride_) {
            *reinterpret_cast<T*>(&mem[i]) = *it++;
        }
        set_data(mem);
        return it;
    }
    /**
     * get data with iterator
     */
    template <typename iterator_type>
    iterator_type get_data(iterator_type const& first) const
    {
        auto data = get_data();
        auto it = first;
        for(size_t i = offset_; i < data.size(); i += stride_) {
            *it++ = *reinterpret_cast<const T*>(&data[i]);
        }
        return it;
    }
    /**
     * set data from lua table
     *
     * @param table lua table containing the data
     */
    virtual void set_lua(luaponte::object table) {
        auto mem = get_memory();
        size_t j = 1;
        for(size_t i = offset_; i < mem.size(); i += stride_) {
            *reinterpret_cast<T*>(&mem[i]) = luaponte::object_cast<T>(table[j++]);
        }
        set_data(mem);
    }
    /**
     * convert data to lua table
     *
     * @param L lua state (passed in by luaponte)
     * @return lua table containing a copy of the data
     */
    virtual luaponte::object get_lua(lua_State *L) {
        auto data = get_data();
        luaponte::object table = luaponte::newtable(L);
        std::size_t j = 1;
        for(size_t i = offset_; i < data.size(); i += stride_) {
            table[j++] = boost::cref(*reinterpret_cast<T*>(&data[i]));
        }
        return table;
    }
protected:
    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     */
    virtual cuda::host::vector<uint8_t> get_memory() const = 0;
    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_data() const = 0;
    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_data(cuda::host::vector<uint8_t> const& memory) = 0;
private:
    /** gpu data stride */
    size_t stride_;
    /** typed data offset */
    size_t offset_;
};

/** typed gpu particle data */
template<typename T>
class particle_data_gpu : public particle_data_typed<T>
{
public:
    /**
     * gpu particle data constructor
     *
     * @param size size of the underlying cuda::vector
     * @param update_function optional update function
     */
    particle_data_gpu(unsigned int size, std::function<void()> update_function)
    : particle_data_typed<T>(sizeof(T), 0), data_(size), update_function_(update_function) {
        if (!update_function_) {
            update_function_ = [](){};
        }
    }

    /**
     * query gpu flag
     *
     * @return true
     */
    virtual bool gpu() const {
        return true;
    }

    /**
     * obtain non-const reference to the stored data
     *
     * @return non-const reference to the cached cuda vector
     *
     * note that this function intentionally does not call
     * the update function, the reason being that the update
     * function will most likely want to use this method
     * to update the data, so calling it here could result
     * in infinite loops
     */
    cache<cuda::vector<T>>& mutable_data() {
        return data_;
    }
    /**
     * obtain const reference to the stored data
     *
     * @return const reference to the cached cuda vector
     */
    cache<cuda::vector<T>> const& data() const {
        update_function_();
        return data_;
    }
protected:
    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     */
    virtual cuda::host::vector<uint8_t> get_memory() const {
        return cuda::host::vector<uint8_t>(data_->size() * sizeof(T));
    }
    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_data() const {
        cuda::host::vector<uint8_t> mem(data_->size() * sizeof(T));
        auto const& g_input = read_cache(data());
        cuda::copy(g_input.begin(), g_input.end(), reinterpret_cast<T*>(&*mem.begin()));
        return mem;
    }
    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_data(cuda::host::vector<uint8_t> const& mem) {
        auto output = make_cache_mutable(mutable_data ());
        auto ptr = reinterpret_cast<T const*>(&*mem.begin());
#ifdef USE_VERLET_DSFUN
        cuda::memset(output->begin(), output->begin() + output->capacity(), 0);
#endif
        cuda::copy(ptr, ptr+output->size(), output->begin());
    }
private:
    /** cached cuda::vector */
    cache<cuda::vector<T>> data_;
    /** optional update function */
    std::function<void()> update_function_;
    /** the host wrapper needs access to get_memory, get_data and set_data */
    template<typename host_type, typename gpu_type, bool is_tuple>
    friend class particle_data_host_wrapper;
};

/**
 * wrapper to access gpu particle data as host vector
 * (e.g. for gpu particle data of type float4 this wrapper
 * can give access to fixed_vector<3,float> values consiting
 * of the first three components of the float4)
 */
template<typename host_type, typename gpu_type, bool is_tuple>
class particle_data_host_wrapper : public particle_data_typed<host_type>
{
public:
    /**
     * host particle data wrapper constructor
     *
     * @param parent gpu particle data to be wrapped
     */
    particle_data_host_wrapper(const std::shared_ptr<particle_data_gpu<gpu_type>> &parent, size_t offset = 0)
            : particle_data_typed<host_type>(sizeof(gpu_type), offset), parent_ (parent) {
    }
    /**
     * query gpu flag
     *
     * @return false
     */
    virtual bool gpu() const {
        return false;
    }
protected:
    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     *
     * depending on the is_tuple template parameter this obtains uninitialized memory from the underlying gpu data
     * or pre-initialized data containing the current contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_memory() const {
        return is_tuple ? parent_->get_data() : parent_->get_memory();
    }
    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_data() const {
        return parent_->get_data();
    }
    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_data(cuda::host::vector<uint8_t> const& memory) {
        return parent_->set_data(memory);
    }
private:
    /** parent gpu data */
    std::shared_ptr<particle_data_gpu<gpu_type>> parent_;
};

// implementations of static members of particle_data

template<typename T>
inline std::shared_ptr<particle_data_gpu<T>> particle_data::create(unsigned int size, std::function<void()> update_function)
{
    return std::make_shared<particle_data_gpu<T>>(size, update_function);
}

template<typename host_type, typename gpu_type>
inline typename std::enable_if<!std::is_same<host_type, gpu_type>::value, std::shared_ptr<particle_data_host_wrapper<host_type, gpu_type>>>::type
particle_data::create_host_wrapper(std::shared_ptr<particle_data_gpu<gpu_type>> const& parent)
{
    return std::make_shared<particle_data_host_wrapper<host_type, gpu_type>>(parent);
}

template<typename host_type, typename gpu_type>
inline typename std::enable_if<std::is_same<host_type, gpu_type>::value, std::shared_ptr<particle_data_gpu<gpu_type>>>::type
particle_data::create_host_wrapper(std::shared_ptr<particle_data_gpu<gpu_type>> const& parent)
{
    return parent;
}

template<typename tuple_type, int field, typename gpu_type>
inline std::shared_ptr<particle_data_host_wrapper<typename std::tuple_element<field, tuple_type>::type, gpu_type, true>>
particle_data::create_packed_wrapper(std::shared_ptr<particle_data_gpu<gpu_type>> const& parent)
{
    typedef typename std::tuple_element<field, tuple_type>::type host_type;
    static constexpr size_t offset = (field == 0) ? 0 : (sizeof(gpu_type) - sizeof(host_type));
    static_assert(std::tuple_size<tuple_type>::value == 2, "invalid tuple");
    return std::make_shared<particle_data_host_wrapper<host_type, gpu_type, true>>(parent, offset);
}

template<typename T>
inline std::shared_ptr<particle_data_typed<T>> particle_data::cast(std::shared_ptr<particle_data> const& ptr)
{
    if(ptr->type() != typeid(T)) {
        throw std::invalid_argument("invalid cast of particle data from " + std::string(ptr->type().name())
                                    + " to " + std::string(typeid(T).name()));
    }
    return std::static_pointer_cast<particle_data_typed<T>>(ptr);
}

template<typename T>
inline std::shared_ptr<particle_data_gpu<T>> particle_data::cast_gpu(std::shared_ptr<particle_data> const& ptr)
{
    if(!ptr->gpu() || ptr->type() != typeid(T)) {
        throw std::invalid_argument("invalid cast of particle data from " + std::string(ptr->type().name())
                                    + " to " + std::string(typeid(T).name()));
    }
    return std::static_pointer_cast<particle_data_gpu<T>>(ptr);
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_DATA_HPP */
