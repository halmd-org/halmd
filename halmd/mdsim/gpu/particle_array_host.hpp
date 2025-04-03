/*
 * Copyright © 2017 Daniel Kirchner
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

#ifndef HALMD_MDSIM_GPU_PARTICLE_ARRAY_HOST_HPP
#define HALMD_MDSIM_GPU_PARTICLE_ARRAY_HOST_HPP

#include <typeinfo>
#include <halmd/io/logger.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/mdsim/gpu/particle_array_gpu.hpp>
#include <halmd/mdsim/gpu/particle_kernel.hpp>
#include <halmd/mdsim/force_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

// stress tensor wrapper
template<typename T>
class stress_tensor_wrapper
        : public T
{
public:
    template<typename... Args>
    stress_tensor_wrapper(Args&&... args) : T(std::forward<Args>(args)...) {}
};

/** host particle array base class */
class particle_array_host_base
{
public:
    virtual ~particle_array_host_base()
    {
    }

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
    virtual luaponte::object get_lua(lua_State* L) const = 0;

    /**
     * get data offset
     */
    virtual size_t offset() const = 0;

    /**
     * get data stride
     */
    virtual size_t stride() const = 0;

    virtual bool coalesced() const = 0;

    virtual std::shared_ptr<particle_array_gpu_base> parent() const = 0;

};

namespace detail {

template<typename T>
struct particle_array_host_helper
{
    typedef T type;
    static T const& get(cuda::memory::host::vector<uint8_t> const& memory, size_t offset)
    {
        return *reinterpret_cast<const T*>(&memory[offset]);
    }
    static void set(cuda::memory::host::vector<uint8_t>& memory, size_t offset, T const& value)
    {
        *reinterpret_cast<T*>(&memory[offset]) = value;
    }
};

/**
 * particle data converter - specialization for stress tensor
 */
template<typename T>
struct particle_array_host_helper<stress_tensor_wrapper<T>>
{
    typedef T type;
    static T get(cuda::memory::host::vector<uint8_t> const& memory, size_t offset)
    {
        unsigned int stride = memory.capacity() / (sizeof(typename T::value_type) * T::static_size);
        return read_stress_tensor<T>(reinterpret_cast<typename T::value_type const*>(&memory[offset]), stride);
    }

    static void set(cuda::memory::host::vector<uint8_t>& memory, size_t offset, T const& value)
    {
        throw std::runtime_error("attempt to write read-only data");
    }
};

} /* namespace detail */

/** host particle array */
template<typename T_>
class particle_array_host : public particle_array_host_base
{
    typedef detail::particle_array_host_helper<T_> helper;
    typedef typename helper::type T;
public:
    particle_array_host(std::shared_ptr<particle_array_gpu_base> const& parent, size_t offset, size_t stride, bool coalesced = false);
    virtual ~particle_array_host();
    virtual std::type_info const& type() const
    {
        return typeid(T);
    }

    static std::shared_ptr<particle_array_host> cast(std::shared_ptr<particle_array_host_base> base);

    /**
     * set data with iterator
     */
    template <typename iterator_type>
    iterator_type set_data(iterator_type const& first)
    {
        auto mem = coalesced_ ? parent_->get_host_data() : parent_->get_host_memory();
        auto it = first;
        // TODO: handle sizes & ghost particles
        size_t offset = offset_;
        for (size_t i = 0; i < parent_->nparticle(); i++) {
            helper::set(mem, offset, *it++);
            offset += stride_;
        }
        parent_->set_host_data(mem);
        return it;
    }

    /**
     * get data with iterator
     */
    template <typename iterator_type>
    iterator_type get_data(iterator_type const& first) const
    {
        auto data = parent_->get_host_data();
        auto it = first;
        // TODO: handle sizes & ghost particles
        size_t offset = offset_;
        for(size_t i = 0; i < parent_->nparticle(); i++) {
            *it++ = helper::get(data, offset);
            offset += stride_;
        }
        return it;
    }

    /**
     * set data from lua table
     *
     * @param table lua table containing the data
     */
    virtual void set_lua(luaponte::object table);

    /**
     * convert data to lua table
     *
     * @param L lua state (passed in by luaponte)
     * @return lua table containing a copy of the data
     */
    virtual luaponte::object get_lua(lua_State *L) const;

    virtual std::shared_ptr<particle_array_gpu_base> parent() const
    {
        return parent_;
    }

    virtual size_t offset() const
    {
        return offset_;
    }

    virtual size_t stride() const
    {
        return stride_;
    }

    virtual bool coalesced() const
    {
        return coalesced_;
    }

private:
    std::shared_ptr<particle_array_gpu_base> parent_;
    size_t offset_;
    size_t stride_;
    bool coalesced_;
};

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_ARRAY_HOST_HPP */
