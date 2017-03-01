/*
 * Copyright © 2016      Daniel Kirchner
 * Copyright © 2010-2012 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2012 Peter Colberg
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

#include <halmd/config.hpp>

#include <halmd/algorithm/gpu/iota.hpp>
#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_kernel.hpp>
#include <halmd/mdsim/gpu/velocity.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <luaponte/out_value_policy.hpp>

#include <algorithm>
#include <exception>
#include <iterator>
#include <numeric>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * Allocate microscopic system state.
 *
 * @param particles number of particles per type or species
 */
template <int dimension, typename float_type>
particle<dimension, float_type>::particle(size_type nparticle, unsigned int nspecies)
  : // allocate global device memory
    nparticle_(nparticle)
  , array_size_((nparticle + 128 - 1) & ~(128 - 1)) // round upwards to multiple of 128
  , nspecies_(std::max(nspecies, 1u))
  // FIXME default CUDA kernel execution dimensions
  , dim_(device::validate(cuda::config(array_size_ / 128, 128)))
  , id_(array_size_)
  , reverse_id_(array_size_)
  // enable auxiliary variables by default to allow sampling of initial state
  , force_zero_(true)
  , force_dirty_(true)
  , aux_dirty_(true)
  , aux_enabled_(true)
{
    // prepare initialization values
    struct {
        fixed_vector<float, 3> position;
        unsigned int species;
    } position_init_value = {
      fixed_vector<float, 3> (0.0f), 0
    }, position_ghost_init_value = {
      fixed_vector<float, 3> (0.0f), -1U
    };
    struct {
        fixed_vector<float, 3> velocity;
        float mass;
    } velocity_init_value = {
      fixed_vector<float, 3> (0.0f), 1.0f
    };
    // register particle arrays
    auto gpu_position_array = gpu_data_["position"] = std::make_shared<particle_array_gpu<gpu_position_type>>
            (dim_, nparticle_, array_size_, position_init_value, position_ghost_init_value);
    auto gpu_image_array = gpu_data_["image"] = std::make_shared<particle_array_gpu<gpu_image_type>>
            (dim_, nparticle_, array_size_);
    auto gpu_velocity_array = gpu_data_["velocity"] = std::make_shared<particle_array_gpu<gpu_velocity_type>>
            (dim_, nparticle_, array_size_, velocity_init_value, velocity_init_value);
    auto gpu_force_array = gpu_data_["force"] = std::make_shared<particle_array_gpu<gpu_force_type>>
            (dim_, nparticle_, array_size_, [this]() { this->update_force_(); });
    auto gpu_en_pot_array = gpu_data_["potential_energy"] = std::make_shared<particle_array_gpu<gpu_en_pot_type>>
            (dim_, nparticle_, array_size_, [this]() { this->update_force_(); });
    // TODO: automatically handle the larger array size for example with an explicit specialization for a stress_tensor_wrapper type
    auto gpu_stress_pot_array = gpu_data_["potential_stress_tensor"] = std::make_shared<particle_array_gpu<gpu_stress_pot_type>>
            (dim_, nparticle_, array_size_ * stress_pot_type::static_size, [this]() { this->update_force_(); });

    // register host data wrappers for packed data
    host_data_["position"] = std::make_shared<particle_array_host<position_type>>(gpu_position_array, 0, sizeof(float4), true);
    host_data_["species"] = std::make_shared<particle_array_host<species_type>>(gpu_position_array, sizeof(float) * 3, sizeof(float) * 4, true);
    host_data_["velocity"] = std::make_shared<particle_array_host<velocity_type>>(gpu_velocity_array, 0, sizeof(float4), true);
    host_data_["mass"] = std::make_shared<particle_array_host<mass_type>>(gpu_velocity_array, sizeof(float) * 3, sizeof(float) * 4, true);

    // register host wrappers for other data
    host_data_["force"] = std::make_shared<particle_array_host<force_type>>(gpu_force_array, 0, sizeof(gpu_force_type));
    host_data_["image"] = std::make_shared<particle_array_host<image_type>>(gpu_image_array, 0, sizeof(gpu_image_type));
    host_data_["potential_energy"] = std::make_shared<particle_array_host<en_pot_type>>(gpu_en_pot_array, 0, sizeof(gpu_en_pot_type));
    host_data_["potential_stress_tensor"] = std::make_shared<particle_array_host<stress_pot_type>>(gpu_stress_pot_array, 0, sizeof(gpu_stress_pot_type));

    {
        auto id = make_cache_mutable(id_);
        auto reverse_id = make_cache_mutable(reverse_id_);
        iota(id->begin(), id->end(), 0);
        iota(reverse_id->begin(), reverse_id->end(), 0);
    }

    LOG_DEBUG("number of CUDA execution blocks: " << dim_.blocks_per_grid());
    LOG_DEBUG("number of CUDA execution threads per block: " << dim_.threads_per_block());

    if (typeid(float_type) == typeid(float)) {
        LOG("integrate using single precision");
    }

    try {
        cuda::copy(nparticle_, get_particle_kernel<dimension, float_type>().nbox);
        cuda::copy(nspecies_, get_particle_kernel<dimension, float_type>().ntype);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy particle parameters to device symbols");
        throw;
    }

    LOG("number of particles: " << nparticle_);
    LOG("number of particle species: " << nspecies_);
    LOG_DEBUG("capacity of data arrays: " << array_size_);
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::aux_enable()
{
    LOG_TRACE("enable computation of auxiliary variables");
    aux_enabled_ = true;
}

/**
 * rearrange particles by permutation
 */
template <int dimension, typename float_type>
void particle<dimension, float_type>::rearrange(cuda::vector<unsigned int> const& g_index)
{
    auto g_position = make_cache_mutable(mutable_data<gpu_position_type>("position"));
    auto g_image = make_cache_mutable(mutable_data<gpu_image_type>("image"));
    auto g_velocity = make_cache_mutable(mutable_data<gpu_velocity_type>("velocity"));

    scoped_timer_type timer(runtime_.rearrange);

    position_array_type position(array_size_);
    image_array_type image(array_size_);
    velocity_array_type velocity(array_size_);
    id_array_type id(array_size_);

    cuda::configure(dim_.grid, dim_.block);
    get_particle_kernel<dimension, float_type>().r.bind(*g_position);
    get_particle_kernel<dimension, float_type>().image.bind(*g_image);
    get_particle_kernel<dimension, float_type>().v.bind(*g_velocity);
    get_particle_kernel<dimension, float_type>().id.bind(read_cache(id_));
    get_particle_kernel<dimension, float_type>().rearrange(g_index, position, image, velocity, id, nparticle_);

    position.swap(*g_position);
    image.swap(*g_image);
    velocity.swap(*g_velocity);
    cuda::copy(id.begin(), id.begin() + id.capacity(), make_cache_mutable(id_)->begin());

    auto reverse_id = make_cache_mutable(reverse_id_);

    iota(reverse_id->begin(), reverse_id->begin() + reverse_id->capacity(), 0);
    radix_sort(id.begin(), id.begin() + nparticle_, reverse_id->begin());
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::update_force_(bool with_aux)
{
    on_prepend_force_();          // ask force modules whether force/aux cache is dirty

    if (force_dirty_ || (with_aux && aux_dirty_)) {
        if (with_aux && aux_dirty_) {
            if (!force_dirty_) {
                LOG_WARNING_ONCE("auxiliary variables inactive in prior force computation, use aux_enable()");
            }
            aux_enabled_ = true;  // turn on computation of aux variables
        }
        LOG_TRACE("request force" << std::string(aux_enabled_ ? " and auxiliary variables" : ""));

        force_zero_ = true;       // tell first force module to reset the force
        on_force_();              // compute forces
        force_dirty_ = false;     // mark force cache as clean
        if (aux_enabled_) {
            aux_dirty_ = false;   // aux cache is clean only if requested
        }
        aux_enabled_ = false;     // disable aux variables for next call
    }
    on_append_force_();
}

template<int dimension, typename float_type>
void particle<dimension, float_type>::insert(std::shared_ptr<particle> const &new_particles)
{
    // TODO
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
}

template <typename particle_type>
static luaponte::object wrap_get(particle_type const& particle, lua_State* L, std::string const& name)
{
    return particle.get_host_array(name)->get_lua(L);
}

template <typename particle_type>
static void wrap_set(particle_type const& particle, std::string const& name, luaponte::object object)
{
    particle.get_host_array(name)->set_lua(object);
}

template <typename T>
static bool equal(std::shared_ptr<T const> self, std::shared_ptr<T const> other)
{
    // compare pointers of managed objects. I could not get working
    // owner_equal() or owner_before() with shared_ptr's passed from Lua
    return self == other;
}

template<typename float_type>
struct variant_name;

template<>
struct variant_name<float>
{
    static constexpr const char *name = "float";
};

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template<>
struct variant_name<dsfloat>
{
    static constexpr const char *name = "dsfloat";
};
#endif

template <int dimension, typename float_type>
void particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string class_name = "particle_" + std::to_string(dimension) + "_" + std::string(variant_name<float_type>::name);
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                class_<particle, std::shared_ptr<particle>>(class_name.c_str())
                    .def(constructor<size_type, unsigned int>())
                    .property("nparticle", &particle::nparticle)
                    .property("array_size", &particle::array_size)
                    .property("nspecies", &particle::nspecies)
                    .def("insert", &particle::insert)
                    .def("get", &wrap_get<particle>)
                    .def("set", &wrap_set<particle>)
                    .def("shift_velocity", &shift_velocity<particle>)
                    .def("shift_velocity_group", &shift_velocity_group<particle>)
                    .def("rescale_velocity", &rescale_velocity<particle>)
                    .def("rescale_velocity_group", &rescale_velocity_group<particle>)
                    .def("shift_rescale_velocity", &shift_rescale_velocity<particle>)
                    .def("shift_rescale_velocity_group", &shift_rescale_velocity_group<particle>)
                    .property("dimension", &wrap_dimension<dimension, float_type>)
                    .def("aux_enable", &particle::aux_enable)
                    .def("on_prepend_force", &particle::on_prepend_force)
                    .def("on_force", &particle::on_force)
                    .def("on_append_force", &particle::on_append_force)
                    .def("__eq", &equal<particle>) // operator= in Lua
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("rearrange", &runtime::rearrange)
                    ]
                    .def_readonly("runtime", &particle::runtime_)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_particle(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    particle<3, float>::luaopen(L);
    particle<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    particle<3, dsfloat>::luaopen(L);
    particle<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class particle<3, float>;
template class particle<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class particle<3, dsfloat>;
template class particle<2, dsfloat>;
#endif

} // namespace gpu
} // namespace mdsim
} // namespace halmd
