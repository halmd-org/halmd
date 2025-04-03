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

#include <halmd/algorithm/gpu/iota.hpp>
#include <halmd/mdsim/gpu/particle_array_gpu.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template<typename T>
particle_array_gpu<T>::~particle_array_gpu()
{}

template<typename T>
struct ValueTypeTrait;

template<>
struct ValueTypeTrait<float> : std::integral_constant<ValueType, ValueType::FLOAT> {};
template<>
struct ValueTypeTrait<float2> : std::integral_constant<ValueType, ValueType::FLOAT2> {};
template<>
struct ValueTypeTrait<float4> : std::integral_constant<ValueType, ValueType::FLOAT4> {};
template<>
struct ValueTypeTrait<dsfloat> : std::integral_constant<ValueType, ValueType::DSFLOAT> {};
template<>
struct ValueTypeTrait<fixed_vector<dsfloat, 2>> : std::integral_constant<ValueType, ValueType::DSFLOAT2> {};
template<>
struct ValueTypeTrait<fixed_vector<dsfloat, 4>> : std::integral_constant<ValueType, ValueType::DSFLOAT4> {};
template<>
struct ValueTypeTrait<unsigned int> : std::integral_constant<ValueType, ValueType::UINT> {};
template<>
struct ValueTypeTrait<uint2> : std::integral_constant<ValueType, ValueType::UINT2> {};
template<>
struct ValueTypeTrait<uint4> : std::integral_constant<ValueType, ValueType::UINT4> {};
template<>
struct ValueTypeTrait<int> : std::integral_constant<ValueType, ValueType::INT> {};
template<>
struct ValueTypeTrait<int2> : std::integral_constant<ValueType, ValueType::INT2> {};
template<>
struct ValueTypeTrait<int4> : std::integral_constant<ValueType, ValueType::INT4> {};

template<typename T>
ValueType particle_array_gpu<T>::value_type() const
{
    return ValueTypeTrait<T>::value;
}

template<typename T>
std::shared_ptr<particle_array_gpu<T>> particle_array_gpu<T>::cast(std::shared_ptr<particle_array_gpu_base> base)
{
    if(base->value_type() != ValueTypeTrait<T>::value) {
        throw std::runtime_error("invalid cast");
    }
    return std::static_pointer_cast<particle_array_gpu<T>>(base);
}

cuda::config get_default_config(size_t n) {
    cuda::device::properties prop(device::get());
    size_t block_size = 128;
    size_t grid_size = n / block_size;
    while (grid_size > prop.max_grid_size().x && block_size <= prop.max_threads_per_block()/2) {
        block_size <<= 1;
        grid_size = (grid_size + 1) >> 1;
    }
    if(grid_size * block_size != n) {
        throw std::runtime_error("misaligned particle array");
    }
    return device::validate(cuda::config(grid_size, block_size));
}

template<typename T>
struct particle_array_gpu_helper
{
    static cuda::memory::host::vector<uint8_t> get_host_memory(
      cache<typename particle_array_gpu<T>::gpu_vector_type> const& data
    )
    {
        cuda::memory::host::vector<uint8_t> mem(data->size() * sizeof(T));
        return mem;
    }

    static cuda::memory::host::vector<uint8_t> get_host_data(
      cache<typename particle_array_gpu<T>::gpu_vector_type> const& data
    )
    {
        auto const& g_input = read_cache(data);
        cuda::memory::host::vector<uint8_t> mem(g_input.size() * sizeof(T));
        mem.reserve(g_input.capacity() * sizeof(T));
        cuda::copy(g_input.begin(), g_input.begin() + g_input.capacity(), reinterpret_cast<T*>(&*mem.begin()));
        return mem;
    }

    static void set_host_data(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , cuda::memory::host::vector<uint8_t> const& mem
    )
    {
        auto output = make_cache_mutable(data);
        auto ptr = reinterpret_cast<T const*>(&*mem.begin());
        cuda::copy(ptr, ptr + (mem.size() / sizeof (T)), output->begin());
    }

    static void initialize_value(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , T const& init_value
      , T const& ghost_init_value
      , unsigned int nparticle
    )
    {
        auto output = make_cache_mutable(data);
        configure_kernel(particle_initialize_wrapper<T>::kernel.initialize, get_default_config(output->size()), true);
        particle_initialize_wrapper<T>::kernel.initialize (output->data(), init_value, ghost_init_value, nparticle);
    }

    static void initialize_zero(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , unsigned int nparticle
    )
    {
        auto output = make_cache_mutable(data);
        cuda::memset(output->begin(), output->end(), 0);
    }

    template<typename U = T>
    static typename std::enable_if<std::is_same<U, unsigned int>::value>::type
    initialize_iota(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , unsigned int nparticle
    )
    {
        auto output = make_cache_mutable(data);
        iota(output->begin(), output->end(), 0);
    }

    template<typename U = T>
    static typename std::enable_if<!std::is_same<U, unsigned int>::value>::type
    initialize_iota(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , unsigned int nparticle
    )
    {
        throw std::runtime_error("iota initialization only supported for integer types");
    }
};

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template<size_t dimension>
struct particle_array_gpu_helper<fixed_vector<dsfloat, dimension>>
{
    typedef typename type_traits<dimension, dsfloat>::gpu::vector_type T;
    typedef typename particle_array_gpu<T>::base_value_type base_value_type;

    static cuda::memory::host::vector<uint8_t> get_host_memory(
        cache<typename particle_array_gpu<T>::gpu_vector_type> const& data
    )
    {
        cuda::memory::host::vector<uint8_t> mem(data->size() * sizeof(base_value_type));
        mem.reserve(data->size() * 2);
        return mem;
    }

    static cuda::memory::host::vector<uint8_t> get_host_data(
        cache<typename particle_array_gpu<T>::gpu_vector_type> const& data
    )
    {
        cuda::memory::device::vector<base_value_type> const& g_input = read_cache(data);
        cuda::memory::host::vector<uint8_t> mem(g_input.size() * sizeof(base_value_type));
        mem.reserve(g_input.capacity() * sizeof(base_value_type));
        cuda::copy(g_input.begin(), g_input.begin() + g_input.capacity(), reinterpret_cast<base_value_type*>(&*mem.begin()));
        return mem;
    }

    static void set_host_data(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data, cuda::memory::host::vector<uint8_t> const& mem
    )
    {
        cuda::memory::device::vector<base_value_type> &output = *make_cache_mutable(data);
        auto ptr = reinterpret_cast<base_value_type const*>(&*mem.begin());
        cuda::copy(ptr, ptr + (mem.capacity() / sizeof (base_value_type)), output.begin());
    }

    static void initialize_value(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , base_value_type const& init_value
      , base_value_type const& ghost_init_value
      , unsigned int nparticle
    )
    {
        auto output = make_cache_mutable(data);
        configure_kernel(dsfloat_particle_initialize_wrapper<dimension>::kernel.initialize, get_default_config(output->size()), true);
        dsfloat_particle_initialize_wrapper<dimension>::kernel.initialize (output->data(), init_value, ghost_init_value, nparticle);
    }

    static void initialize_zero(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , unsigned int nparticle
    )
    {
        cuda::memory::device::vector<base_value_type> &output = *make_cache_mutable(data);
        cuda::memset(output.begin(), output.begin() + output.capacity(), 0);
    }

    static void initialize_iota(
        cache<typename particle_array_gpu<T>::gpu_vector_type>& data
      , unsigned int nparticle
    )
    {
        throw std::runtime_error("iota initialization only supported for integer types");
    }
};

template<>
struct particle_array_gpu_helper<dsfloat> : particle_array_gpu_helper<fixed_vector<dsfloat, 1>> {};
#endif // USE_GPU_DOUBLE_SINGLE_PRECISION

template<typename T>
cuda::memory::host::vector<uint8_t> particle_array_gpu<T>::get_host_memory() const
{
    return particle_array_gpu_helper<T>::get_host_memory(data_);
}

template<typename T>
cuda::memory::host::vector<uint8_t> particle_array_gpu<T>::get_host_data() const
{
    update_function_();
    return particle_array_gpu_helper<T>::get_host_data(data_);
}

template<typename T>
void particle_array_gpu<T>::set_host_data(cuda::memory::host::vector<uint8_t> const& mem)
{
    return particle_array_gpu_helper<T>::set_host_data(data_, mem);
}

template<typename T>
void particle_array_gpu<T>::initialize()
{
    switch(init_type_) {
        case InitType::ZERO:
            particle_array_gpu_helper<T>::initialize_zero(data_, nparticle_);
            break;
        case InitType::IOTA:
            particle_array_gpu_helper<T>::initialize_iota(data_, nparticle_);
            break;
        case InitType::VALUE:
            particle_array_gpu_helper<T>::initialize_value(data_, init_value_, ghost_init_value_, nparticle_);
            break;
    }
}

template class particle_array_gpu<float>;
template class particle_array_gpu<float2>;
template class particle_array_gpu<float4>;
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class particle_array_gpu<dsfloat>;
template class particle_array_gpu<fixed_vector<dsfloat, 2>>;
template class particle_array_gpu<fixed_vector<dsfloat, 4>>;
#endif
template class particle_array_gpu<unsigned int>;
template class particle_array_gpu<uint2>;
template class particle_array_gpu<uint4>;
template class particle_array_gpu<int>;
template class particle_array_gpu<int2>;
template class particle_array_gpu<int4>;

} // namespace gpu
} // namespace mdsim
} // namespace halmd
