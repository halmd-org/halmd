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

#ifndef HALMD_MDSIM_GPU_PARTICLE_ARRAY_GPU_HPP
#define HALMD_MDSIM_GPU_PARTICLE_ARRAY_GPU_HPP

#include <typeinfo>
#include <halmd/io/logger.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/gpu/dsfloat_cuda_vector.hpp>
#include <halmd/mdsim/gpu/particle_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

enum class ValueType : unsigned int {
    FLOAT,
    FLOAT2,
    FLOAT4,
    DSFLOAT,
    DSFLOAT2,
    DSFLOAT4,
    UINT,
    UINT2,
    UINT4,
    INT,
    INT2,
    INT4
};

enum class InitType : unsigned int {
    VALUE,
    ZERO,
    IOTA
};

/** gpu particle array base class */
class particle_array_gpu_base
{
public:
    virtual ~particle_array_gpu_base()
    {
    }
    virtual ValueType value_type() const = 0;
    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     */
    virtual cuda::memory::host::vector<uint8_t> get_host_memory() const = 0;

    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::memory::host::vector<uint8_t> get_host_data() const = 0;

    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_host_data(cuda::memory::host::vector<uint8_t> const& memory) = 0;

    /**
     * query cache observer
     *
     * @return a cache observer reflecting the current state of the cache
     */
    virtual cache<> cache_observer() const = 0;

    /**
     * return number of particles
     */
    virtual size_t nparticle() const = 0;
};

template<typename T>
struct particle_array_gpu_traits
{
    typedef T base_value_type;
    typedef cuda::memory::device::vector<T> gpu_vector_type;
};

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template<size_t dimension>
struct particle_array_gpu_traits<fixed_vector<dsfloat, dimension> >
{
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type base_value_type;
    typedef dsfloat_vector<cuda::memory::device::vector<base_value_type>> gpu_vector_type;
};

template<>
struct particle_array_gpu_traits<dsfloat> : particle_array_gpu_traits<fixed_vector<dsfloat, 1>> {};
#endif // USE_GPU_DOUBLE_SINGLE_PRECISION


template<typename T>
class particle_array_gpu : public particle_array_gpu_base
{
public:
    typedef typename particle_array_gpu_traits<T>::base_value_type base_value_type;
    typedef typename particle_array_gpu_traits<T>::gpu_vector_type gpu_vector_type;

    template<typename init_type, typename ghost_init_type>
    particle_array_gpu(
        cuda::config const& dim
      , unsigned int nparticle
      , unsigned int size
      , init_type const& init_value
      , ghost_init_type const& ghost_init_value
      , std::function<void()> update_function = std::function<void()>()
    ) : data_(size)
      , update_function_(update_function)
      , init_type_(InitType::VALUE)
      , nparticle_(nparticle)
    {
        // due to a bug in clang moving this to the default argument initialization
        // above results in a compiler crash
        if (!update_function_) {
            update_function_ = []{};
        }
        static_assert(sizeof(init_type) == sizeof(base_value_type), "invalid size of initialization value");
        static_assert(sizeof(ghost_init_type) == sizeof(base_value_type), "invalid size of ghost initialization value");
        memcpy(&init_value_, &init_value, sizeof(base_value_type));
        memcpy(&ghost_init_value_, &ghost_init_value, sizeof(base_value_type));
        initialize();
    }

    particle_array_gpu(
        cuda::config const& dim
      , unsigned int nparticle
      , unsigned int size
      , std::function<void()> update_function = std::function<void()>()
    ) : data_(size)
      , update_function_(update_function)
      , init_type_(InitType::ZERO)
      , nparticle_(nparticle)
    {
        // due to a bug in clang moving this to the default argument initialization
        // above results in a compiler crash
        if (!update_function_) {
            update_function_ = []{};
        }
        initialize();
    }

    static std::shared_ptr<particle_array_gpu> cast(std::shared_ptr<particle_array_gpu_base> base);

    virtual ~particle_array_gpu();
    virtual ValueType value_type() const;

    cache<gpu_vector_type>& mutable_data()
    {
        return data_;
    }

    cache<gpu_vector_type> const& data() const
    {
        update_function_();
        return data_;
    }

    /**
     * get memory
     *
     * @return a page-locked memory vector to be filled with data and passed to set_data
     */
    virtual cuda::memory::host::vector<uint8_t> get_host_memory() const;

    /**
     * get data
     *
     * @return a page-locked memory vector containing the contents of the underlying gpu data
     */
    virtual cuda::memory::host::vector<uint8_t> get_host_data() const;

    /**
     * set data
     *
     * @param memory page-locked memory vector containing data to be copied to the underlying gpu data
     *               should have been obtained with get_memory
     */
    virtual void set_host_data(cuda::memory::host::vector<uint8_t> const& mem);

    /**
     * return number of particles
     */
    virtual size_t nparticle() const
    {
        return nparticle_;
    }

    /**
     * query cache observer
     *
     * @return a cache observer reflecting the current state of the cache
     */
    virtual cache<> cache_observer() const
    {
        return data_;
    }

private:
    void initialize();

    cache<gpu_vector_type> data_;
    std::function<void()> update_function_;
    InitType init_type_;
    base_value_type init_value_;
    base_value_type ghost_init_value_;
    size_t nparticle_;
};

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_ARRAY_GPU_HPP */
