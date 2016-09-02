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

#ifndef HALMD_MDSIM_GPU_PARTICLE_ARRAY_HPP
#define HALMD_MDSIM_GPU_PARTICLE_ARRAY_HPP

#include <typeinfo>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

/*
 * wrapper templates used to mark special types that need a custom particle_array_typed implementation
 * the wrapper template is necessary to create template specializations
 */
// stress tensor wrapper
template<typename T>
class stress_tensor_wrapper : public T {
public:
    template<typename... Args>
    stress_tensor_wrapper(Args&&... args) : T(std::forward<Args>(args)...) {
    }
};

// forward declarations
template<typename T>
class particle_array_typed;
template<typename T>
class particle_array_gpu;
template<typename host_type>
class particle_array_host_wrapper;

/** particle array base class */
class particle_array : public std::enable_shared_from_this<particle_array>
{
public:
    /**
     * create gpu particle array of given type
     *
     * @param size number of particles
     * @param update_function optional update function
     * @return shared pointer to the new particle array
     */
    template<typename T>
    static inline std::shared_ptr<particle_array_gpu<T>> create(unsigned int size, std::function<void()> update_function = std::function<void()>());

    /**
     * cast generic particle array to typed data
     *
     * @param ptr generic particle array
     * @return typed particle array
     *
     * throws an exception if the type of the array does not match
     */
    template<typename T>
    static inline std::shared_ptr<particle_array_typed<T>> cast(std::shared_ptr<particle_array> const& ptr);

    /**
     * cast generic particle array to typed gpu array
     *
     * @param ptr generic particle array
     * @return typed gpu particle array
     *
     * throws an exception if the type of the array does not match
     */
    template<typename T>
    static inline std::shared_ptr<particle_array_gpu<T>> cast_gpu(std::shared_ptr<particle_array> const& ptr);

    /**
     * cast generic particle array to typed host wrapper array
     *
     * @param ptr generic particle array
     * @return typed host wrapper particle array
     *
     * throws an exception if the type of the array does not match
     */
    template<typename T>
    static inline std::shared_ptr<particle_array_host_wrapper<T>> cast_host_wrapper(std::shared_ptr<particle_array> const& ptr);

    /**
     * create host wrapper for particle array of given gpu and host type
     *
     * @param parent gpu particle array
     * @return shared pointer to the particle array wrapper
     *
     * this version returns an actual wrapper object for the case that gpu_type and host_type
     * are actually different
     */
    template<typename host_type, typename gpu_type>
    static inline typename std::enable_if<!std::is_same<host_type, gpu_type>::value,
            std::shared_ptr<particle_array_host_wrapper<host_type>>>::type
    create_host_wrapper(std::shared_ptr<particle_array_gpu<gpu_type>> const& parent);

    /**
     * create host wrapper for particle array of given gpu and host type
     *
     * @param parent gpu particle array
     * @return shared pointer to the particle array wrapper
     *
     * this version directly passes down the parent object for the case that gpu_type and host_type
     * are actually the same (e.g. in the case of scalar values)
     */
    template<typename host_type, typename gpu_type>
    static inline typename std::enable_if<std::is_same<host_type, gpu_type>::value,
            std::shared_ptr<particle_array_gpu<gpu_type>>>::type
    create_host_wrapper(std::shared_ptr<particle_array_gpu<gpu_type>> const& parent);

    /**
     * create host wrapper for an array of tuple elements in a packed gpu particle array
    *
    * @param parent gpu particle array
    * @return shared pointer to the particle array wrapper
    */
    template<typename tuple_type, int field, typename gpu_type>
    static inline std::shared_ptr<particle_array_host_wrapper<typename std::tuple_element<field, tuple_type>::type>>
    create_packed_wrapper(std::shared_ptr<particle_array_gpu<gpu_type>> const& parent);

    /**
     * empty virtual destructor
     *
     * ensures that the particle array is correctly destroyed even if stored
     * as (shared) pointer to the base class
     */
    virtual ~particle_array(){}

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
     * query gpu flag
     *
     * @return true if the particle array contains raw gpu data, false otherwise
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
    virtual luaponte::object get_lua(lua_State* L) const = 0;

    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     */
    virtual cuda::host::vector<uint8_t> get_gpu_memory() const = 0;
    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_gpu_data() const = 0;
    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_gpu_data(cuda::host::vector<uint8_t> const& memory) = 0;
};

namespace detail {

/**
 * particle data converter
 *
 * converts data copied from GPU storage to the desired data type
 *
 * the default implementation just casts and copies the data directly,
 * but custom specializations can perform arbitrary transformations
 * on the data for specific types
 */
template<typename T>
class particle_data_converter {
public:
    T const& get(cuda::host::vector<uint8_t> const& memory, size_t offset) const {
        return *reinterpret_cast<const T*>(&memory[offset]);
    }
    void set(cuda::host::vector<uint8_t>& memory, size_t offset, T const& value) {
        *reinterpret_cast<T*>(&memory[offset]) = value;
    }
};

/**
 * particle data converter - specialization for stress tensor
 */
template<typename T>
class particle_data_converter<stress_tensor_wrapper<T>> {
public:
    T get(cuda::host::vector<uint8_t> const& memory, size_t offset) const {
        unsigned int stride = memory.capacity() / (sizeof(typename T::value_type) * T::static_size);
        return read_stress_tensor<T>(reinterpret_cast<typename T::value_type const*>(&memory[offset]), stride);
    }
    void set(cuda::host::vector<uint8_t>& memory, size_t offset, T const& value) {
        throw std::runtime_error("attempt to write read-only data");
    }
};

} /* namespace detail */

/** particle array base class with type information */
template<typename T>
class particle_array_typed : public particle_array
{
public:
    /**
     * typed particle array constructor
     *
     * @param stride stride of the underlying gpu data
     * @param offset offset of the typed data within the underlying gpu data
     * @param n_elems number of elements in the array
     */
    particle_array_typed(size_t stride, size_t offset, size_t n_elems)
    : stride_(stride), offset_(offset), n_elems_(n_elems) {
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
        auto mem = get_gpu_memory();
        auto it = first;
        for(size_t i = offset_; i < mem.size(); i += stride_) {
            converter.set(mem, i, *it++);
        }
        // TODO: handle dsfloat data correctly here for tuple arrays
        set_gpu_data(mem);
        return it;
    }
    /**
     * get data with iterator
     */
    template <typename iterator_type>
    iterator_type get_data(iterator_type const& first) const
    {
        auto data = get_gpu_data();
        auto it = first;
        for(size_t i = offset_; i < data.size(); i += stride_) {
            *it++ = converter.get(data, i);
        }
        return it;
    }

    /**
     * set data from lua table
     *
     * @param table lua table containing the data
     */
    virtual void set_lua(luaponte::object table) {
        auto mem = get_gpu_memory();
        size_t j = 1;
        for(size_t i = offset_; i < mem.size(); i += stride_) {
            converter.set(mem, i, luaponte::object_cast<T>(table[j++]));
        }
        set_gpu_data(mem);
    }

    /**
     * convert data to lua table
     *
     * @param L lua state (passed in by luaponte)
     * @return lua table containing a copy of the data
     */
    virtual luaponte::object get_lua(lua_State *L) const {
        auto data = get_gpu_data();
        luaponte::object table = luaponte::newtable(L);
        std::size_t j = 1;
        for(size_t i = offset_; i < data.size(); i += stride_) {
            auto&& value = converter.get(data, i);
            table[j++] = boost::cref(value);
        }
        return table;
    }

    size_t stride() const {
        return stride_;
    }

    size_t offset() const {
        return offset_;
    }

    size_t n_elems() const {
        return n_elems_;
    }
private:
    /** gpu data stride */
    size_t stride_;
    /** typed data offset */
    size_t offset_;
    /** number of elements */
    size_t n_elems_;
    /** particle data converter */
    detail::particle_data_converter<T> converter;
};

/** typed gpu particle array */
template<typename T>
class particle_array_gpu : public particle_array_typed<T>
{
public:
    /**
     * gpu particle array constructor
     *
     * @param size size of the underlying cuda::vector
     * @param update_function optional update function
     */
    particle_array_gpu(unsigned int size, std::function<void()> update_function)
    : particle_array_typed<T>(sizeof(T), 0, size), data_(size), update_function_(update_function) {
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
    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     */
    virtual cuda::host::vector<uint8_t> get_gpu_memory() const {
        return cuda::host::vector<uint8_t>(data_->size() * sizeof(T));
    }
    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_gpu_data() const {
        auto const& g_input = read_cache(data());
        cuda::host::vector<uint8_t> mem(g_input.size() * sizeof(T));
        mem.reserve(g_input.capacity() * sizeof(T));
        cuda::copy(g_input.begin(), g_input.begin()+g_input.capacity(), reinterpret_cast<T*>(&*mem.begin()));
        return mem;
    }
    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_gpu_data(cuda::host::vector<uint8_t> const& mem) {
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
};

/**
 * wrapper to access gpu particle array as host vector
 * (e.g. for gpu particle data of type float4 this wrapper
 * can give access to fixed_vector<3,float> values consiting
 * of the first three components of the float4)
 */
template<typename host_type>
class particle_array_host_wrapper : public particle_array_typed<host_type>
{
public:
    /**
     * host particle array wrapper constructor
     *
     * @param parent gpu particle array to be wrapped
     */
    template<typename gpu_type>
    particle_array_host_wrapper(const std::shared_ptr<particle_array_gpu<gpu_type>> &parent, size_t offset, bool is_tuple)
            : particle_array_typed<host_type>(sizeof(gpu_type), offset, parent->data()->size()), is_tuple_(is_tuple), parent_(parent) {
    }
    /**
     * query gpu flag
     *
     * @return false
     */
    virtual bool gpu() const {
        return false;
    }

    /**
    * query cache observer
    *
    * @return a cache observer reflecting the current state of the cache
    */
    virtual cache<> cache_observer() const {
        return parent_->cache_observer();
    }

    std::shared_ptr<particle_array> parent() {
        return parent_;
    }

    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     *
     * depending on the is_tuple template parameter this obtains uninitialized memory from the underlying gpu data
     * or pre-initialized data containing the current contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_gpu_memory() const {
        return is_tuple_ ? parent_->get_gpu_data() : parent_->get_gpu_memory();
    }
    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::host::vector<uint8_t> get_gpu_data() const {
        return parent_->get_gpu_data();
    }
    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_gpu_data(cuda::host::vector<uint8_t> const& memory) {
        return parent_->set_gpu_data(memory);
    }
    bool is_tuple() const {
        return is_tuple_;
    }
private:
    /** tuple flag */
    bool is_tuple_;
    /** parent gpu data */
    std::shared_ptr<particle_array> parent_;
};

// implementations of static members of particle_array

template<typename T>
inline std::shared_ptr<particle_array_gpu<T>> particle_array::create(unsigned int size, std::function<void()> update_function)
{
    return std::make_shared<particle_array_gpu<T>>(size, update_function);
}

template<typename host_type, typename gpu_type>
inline typename std::enable_if<!std::is_same<host_type, gpu_type>::value, std::shared_ptr<particle_array_host_wrapper<host_type>>>::type
particle_array::create_host_wrapper(std::shared_ptr<particle_array_gpu<gpu_type>> const& parent)
{
    return std::make_shared<particle_array_host_wrapper<host_type>>(parent, 0, false);
}

template<typename host_type, typename gpu_type>
inline typename std::enable_if<std::is_same<host_type, gpu_type>::value, std::shared_ptr<particle_array_gpu<gpu_type>>>::type
particle_array::create_host_wrapper(std::shared_ptr<particle_array_gpu<gpu_type>> const& parent)
{
    return parent;
}

template<typename tuple_type, int field, typename gpu_type>
inline std::shared_ptr<particle_array_host_wrapper<typename std::tuple_element<field, tuple_type>::type>>
particle_array::create_packed_wrapper(std::shared_ptr<particle_array_gpu<gpu_type>> const& parent)
{
    typedef typename std::tuple_element<field, tuple_type>::type host_type;
    static constexpr size_t offset = (field == 0) ? 0 : (sizeof(gpu_type) - sizeof(host_type));
    static_assert(std::tuple_size<tuple_type>::value == 2, "invalid tuple");
    return std::make_shared<particle_array_host_wrapper<host_type>>(parent, offset, true);
}

template<typename T>
inline std::shared_ptr<particle_array_typed<T>> particle_array::cast(std::shared_ptr<particle_array> const& ptr)
{
    if(ptr->type() != typeid(T)) {
        throw std::invalid_argument("invalid cast of particle array from " + std::string(ptr->type().name())
                                    + " to " + std::string(typeid(T).name()));
    }
    return std::static_pointer_cast<particle_array_typed<T>>(ptr);
}

template<typename T>
inline std::shared_ptr<particle_array_gpu<T>> particle_array::cast_gpu(std::shared_ptr<particle_array> const& ptr)
{
    if(!ptr->gpu() || ptr->type() != typeid(T)) {
        throw std::invalid_argument("invalid cast of particle array from " + std::string(ptr->type().name())
                                    + " to " + std::string(typeid(T).name()));
    }
    return std::static_pointer_cast<particle_array_gpu<T>>(ptr);
}

template<typename T>
inline std::shared_ptr<particle_array_host_wrapper<T>> particle_array::cast_host_wrapper(std::shared_ptr<particle_array> const& ptr)
{
    if(ptr->gpu() || ptr->type() != typeid(T)) {
        throw std::invalid_argument("invalid cast of particle array from " + std::string(ptr->type().name())
                                    + " to " + std::string(typeid(T).name()));
    }
    return std::static_pointer_cast<particle_array_host_wrapper<T>>(ptr);
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_ARRAY_HPP */
