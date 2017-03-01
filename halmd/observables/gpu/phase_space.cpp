/*
 * Copyright © 2016       Daniel Kirchner
 * Copyright © 2008-2012  Felix Höfling
 * Copyright © 2008-2012  Peter Colberg
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
#include <halmd/observables/gpu/samples/sample.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/**
 * phase space host cache
 *
 * caches the GPU data of a particle instance in page-locked host memory
 *
 * phase_space creates a host cache for every GPU array the first time it
 * is sampled and stores it in an associative map to avoid duplicates in
 * case that multiple host wrapper arrays backed by the same GPU array
 * are sampled
 */
class phase_space_host_cache {
public:
    phase_space_host_cache(std::shared_ptr<mdsim::gpu::particle_array_gpu_base const> array)
      : array_(array) {}

    bool up_to_date() const
    {
        return array_->cache_observer() == cache_observer_;
    }

    cuda::host::vector<uint8_t>& acquire(void)
    {
        if(!(array_->cache_observer() == cache_observer_)) {
            data_ = array_->get_host_data();
            cache_observer_ = array_->cache_observer();
        }
        return data_;
    }

private:
    cuda::host::vector<uint8_t> data_;
    cache<> cache_observer_;
    std::shared_ptr<mdsim::gpu::particle_array_gpu_base const> array_;
};

/**
 * phase space sampler implementation for typed host samples
 *
 * copies the GPU data to a host sample and provides the actual implementation
 * for the sample related interface of phase_space
 */
template<int dimension, typename scalar_type>
class phase_space_sampler_typed
  : public phase_space_sampler_host
{
public:
    typedef host::samples::sample<dimension, scalar_type> sample_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::gpu::particle_array_host<typename sample_type::data_type> particle_array_type;
    typedef std::map<mdsim::gpu::particle_array_gpu_base*, std::shared_ptr<phase_space_host_cache>> host_cache_type;

    /**
     * Creates a sampler for the given particle group and particle array.
     *
     * @param group         particle group passed down from phase_space
     * @param array         particle array containing the data to be sampled
     * @param nthreads      number of GPU threads, currently needed for dsfloat data
     * @param host_cache    reference to the cache of gpu data stored in phase_space
     */
    phase_space_sampler_typed(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<mdsim::gpu::particle_array_host_base> array
      , std::size_t nthreads
      , host_cache_type& host_cache
    )
      : array_(mdsim::gpu::particle_array_host<typename sample_type::data_type>::cast(array))
      , particle_group_(group)
      , nthreads_(nthreads)
    {
        // find host cache for the gpu particle array backing the queried host wrapper array
        // or create and store a new one if none was created so far
        auto parent = array_->parent();
        auto it = host_cache.find(parent.get());
        if (it == host_cache.end()) {
            host_cache_ = host_cache[parent.get()] = std::make_shared<phase_space_host_cache>(parent);
        } else {
            host_cache_ = it->second;
        }
    }

    /**
     * Static wrapper to directly construct a shared_ptr.
     */
    static std::shared_ptr<phase_space_sampler_typed> create(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<mdsim::gpu::particle_array_host_base> array
      , std::size_t nthreads
      , host_cache_type& host_cache
    )
    {
        return std::make_shared<phase_space_sampler_typed>(group, array, nthreads, host_cache);
    }

    /**
     * acquire a sample
     *
     * copies the data from the cached GPU data to a new sample if the data is not up to date,
     * returns the stored sample otherwise
     */
    virtual std::shared_ptr<sample_base> acquire() {
        typedef typename sample_type::data_type data_type;

        if (!host_cache_->up_to_date() || group_observer_ != particle_group_->ordered()) {
            auto const& data = host_cache_->acquire();
            auto const& group = particle_group_->ordered_host_cached();

            sample_ = std::make_shared<sample_type>(group.size());
            auto& sample_data = sample_->data();
            auto offset = array_->offset();
            auto stride = array_->stride();

            size_t id = 0;
            for (size_t i : group) {
                sample_data[id++] = *reinterpret_cast<data_type const*>(&data[offset + i * stride]);
            }

            group_observer_ = particle_group_->ordered();
        }
        return sample_;
    }

    /**
     * copies the data of a sample to the GPU particle array
     */
    virtual void set(std::shared_ptr<sample_base const> sample_)
    {
        typedef typename sample_type::data_type data_type;

        if(sample_->gpu() || sample_->type() != typeid(typename sample_type::data_type)) {
            throw std::runtime_error("invalid sample data type");
        }
        auto sample = std::static_pointer_cast<sample_type const>(sample_);
        auto const& sample_data = sample->data();

        auto &data = host_cache_->acquire(); // TODO: only for coalesced?
        auto const& group = particle_group_->ordered_host_cached();
        auto offset = array_->offset();
        auto stride = array_->stride();

        size_t id = 0;
        for(size_t i : group) {
            *reinterpret_cast<data_type*>(&data[offset + i * stride]) = sample_data[id++];
        }
        // TODO: only needed for dsfloat types
        if(data.capacity() >= data.size() + nthreads_ * stride) {
            offset += nthreads_ * stride;
            for(size_t i : group) {
                *reinterpret_cast<data_type*>(&data[offset + i * stride]) = data_type(0);
            }
        }
        array_->parent()->set_host_data(data);
    }

    /**
     * returns a lua slot function to be used to acquire a host sample
     */
    virtual luaponte::object acquire_lua(lua_State* L, std::shared_ptr<phase_space_sampler_host> self_)
    {
        auto self = std::static_pointer_cast<phase_space_sampler_typed>(self_);
        std::function<std::shared_ptr<sample_type const>()> fn = [self]() -> std::shared_ptr<sample_type const>
        {
            return std::static_pointer_cast<sample_type const>(self->acquire());
        };

        luaponte::default_converter<std::function<std::shared_ptr<sample_type const>()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }

    /**
     * returns a lua slot function to be used to directly acquire the data of a host sample
     */
    virtual luaponte::object data_lua(lua_State* L, std::shared_ptr<phase_space_sampler_host> self)
    {
        std::function<typename sample_type::array_type const&()> fn = [self]() -> typename sample_type::array_type const&
        {
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
    virtual void set_lua(luaponte::object sample)
    {
        set(luaponte::object_cast<std::shared_ptr<sample_type const>>(sample));
    }

protected:
    /** particle array containing the data to be sampled */
    std::shared_ptr<particle_array_type> array_;
    /** particle group containing the data indices to be sampled */
    std::shared_ptr<particle_group_type> particle_group_;
    /** host cache of the GPU data backing the sample data */
    std::shared_ptr<phase_space_host_cache> host_cache_;
    /** cached sample */
    std::shared_ptr<sample_type> sample_;
    /** cache observer for the index list */
    cache<> group_observer_;
    /** number of GPU threads (currently only used for dsfloat data) */
    std::size_t nthreads_;
};

/**
 * associative map from typeid's to typed phase space sampler create functions
 *
 * used to create the correct phase space sampler based on the typeid of a particle array
 */
static const std::unordered_map<
    std::type_index
  , std::function<std::shared_ptr<phase_space_sampler_host>(
        std::shared_ptr<mdsim::gpu::particle_group>
      , std::shared_ptr<mdsim::gpu::particle_array_host_base>
      , std::size_t
      , std::map<mdsim::gpu::particle_array_gpu_base*, std::shared_ptr<phase_space_host_cache>>&
    )>
> phase_space_sampler_typed_create_map = {
    { typeid(float), phase_space_sampler_typed<1, float>::create }
  , { typeid(fixed_vector<float, 2>), phase_space_sampler_typed<2, float>::create }
  , { typeid(fixed_vector<float, 3>), phase_space_sampler_typed<3, float>::create }
  , { typeid(fixed_vector<float, 4>), phase_space_sampler_typed<4, float>::create }
/*
  , { typeid(double), phase_space_sampler_typed<1, double>::create }
  , { typeid(fixed_vector<double, 2>), phase_space_sampler_typed<2, double>::create }
  , { typeid(fixed_vector<double, 3>), phase_space_sampler_typed<3, double>::create }
  , { typeid(fixed_vector<double, 4>), phase_space_sampler_typed<4, double>::create }
*/
  , { typeid(int), phase_space_sampler_typed<1, int>::create }
  , { typeid(fixed_vector<int, 2>), phase_space_sampler_typed<2, int>::create }
  , { typeid(fixed_vector<int, 3>), phase_space_sampler_typed<3, int>::create }
  , { typeid(fixed_vector<int, 4>), phase_space_sampler_typed<4, int>::create }

  , { typeid(unsigned int), phase_space_sampler_typed<1, unsigned int>::create }
  , { typeid(fixed_vector<unsigned int, 2>), phase_space_sampler_typed<2, unsigned int>::create }
  , { typeid(fixed_vector<unsigned int, 3>), phase_space_sampler_typed<3, unsigned int>::create }
  , { typeid(fixed_vector<unsigned int, 4>), phase_space_sampler_typed<4, unsigned int>::create }
};


/**
 * specialized phase_space sampler for host position data
 *
 * does the same as the generic sampler, but additionally reduces/extends periodic positions
 */
template<int dimension, typename float_type, typename scalar_type>
class phase_space_sampler_position
  : public phase_space_sampler_typed<dimension, scalar_type>
{
public:
    typedef host::samples::sample<dimension, scalar_type> sample_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle_array_host<typename sample_type::data_type> particle_array_type;

    /**
     * Creates a position sampler for the given particle group.
     *
     * @param group             particle group passed down from phase_space
     * @param box               box for periodic border conditions
     * @param position_array    particle array containing the position data
     * @param image_array       particle array containing the image data
     * @param dim               CUDA configuration used to launch CUDA kernels
     * @param host_cache        reference to the cache of gpu data stored in phase_space
     */
    phase_space_sampler_position(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<mdsim::gpu::particle_array_host_base> position_array
      , std::shared_ptr<mdsim::gpu::particle_array_host_base> image_array
      , cuda::config const& dim
      , std::map<mdsim::gpu::particle_array_gpu_base*, std::shared_ptr<phase_space_host_cache>>& host_cache
    )
      : phase_space_sampler_typed<dimension, scalar_type>(group, position_array, dim.threads(), host_cache)
      , box_(box)
      , image_array_(mdsim::gpu::particle_array_host<typename sample_type::data_type>::cast(image_array))
      , dim_(dim)
    {
        auto image_parent = image_array_->parent();
        auto it = host_cache.find(image_parent.get());
        if (it == host_cache.end()) {
            image_host_cache_ = host_cache[image_parent.get()] = std::make_shared<phase_space_host_cache>(image_parent);
        } else {
            image_host_cache_ = it->second;
        }
    }

    /**
     * Static wrapper to directly construct a shared_ptr.
     */
    static std::shared_ptr<phase_space_sampler_position> create(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<mdsim::gpu::particle_array_host_base> position_array
      , std::shared_ptr<mdsim::gpu::particle_array_host_base> image_array
      , cuda::config const& dim
      , std::map<mdsim::gpu::particle_array_gpu_base*, std::shared_ptr<phase_space_host_cache>>& host_cache)
    {
        return std::make_shared<phase_space_sampler_position>(group, box, position_array, image_array, dim, host_cache);
    }

    /**
     * acquire a sample
     *
     * Copies the data from the cached GPU data to a new sample if the data is not up to date,
     * returns the stored sample otherwise.
     * Periodically extends position data.
     */
    virtual std::shared_ptr<sample_base> acquire()
    {
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

            size_t id = 0;
            for (size_t i : group) {
                auto& r = sample_position[id++];
                r = *reinterpret_cast<typename sample_type::data_type const*>(&particle_position[position_offset + i * position_stride]);
                box_->extend_periodic(r, *reinterpret_cast<typename sample_type::data_type const*>(&particle_image[image_offset + i * image_stride]));
            }
        }
        return this->sample_;
    }

    /**
     * Copies the data of a sample to the GPU particle array.
     * This is done by a CUDA kernel that automatically periodically reduces
     * the position data.
     */
    virtual void set(std::shared_ptr<sample_base const> sample)
    {
        // set position data the same as any other kind of data
        phase_space_sampler_typed<dimension, scalar_type>::set(sample);

        // reduce positions on GPU
        typedef typename mdsim::type_traits<dimension, scalar_type>::gpu::coalesced_vector_type gpu_vector_type;
        auto position = make_cache_mutable(mdsim::gpu::particle_array_gpu<
          typename mdsim::gpu::particle<dimension, float_type>::gpu_hp_vector_type
          >::cast(this->array_->parent())->mutable_data());
        auto image = make_cache_mutable(mdsim::gpu::particle_array_gpu<gpu_vector_type>::cast(this->image_array_->parent())->mutable_data());
        auto const& group = read_cache(this->particle_group_->ordered());
        try {
            phase_space_wrapper<dimension>::kernel.r.bind(*position);
            cuda::configure(dim_.grid, dim_.block);
            phase_space_wrapper<dimension>::kernel.reduce_periodic(
                &*group.begin()
              , position->data() // TODO: is this correct for dsfloats?
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
    /** box for periodic border conditions */
    std::shared_ptr<box_type const> box_;
    /** particle array storing the image data */
    std::shared_ptr<particle_array_type> image_array_;
    /** host cache of the GPU data backing the sample data */
    std::shared_ptr<phase_space_host_cache> image_host_cache_;
    /** CUDA configuration for launching kernels */
    cuda::config const dim_;
};

/**
 * phase space sampler implementation for typed gpu samples
 *
 * copies the GPU data to a GPU sample and provides the actual implementation
 * for the sample related interface of phase_space
 */
template<int dimension, typename input_data_type, typename sample_data_type = input_data_type>
class phase_space_sampler_gpu_typed
  : public phase_space_sampler_gpu
{
    typedef phase_space_sample_wrapper<input_data_type, sample_data_type> sample_wrapper;
public:
    typedef samples::sample<dimension, sample_data_type> sample_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::gpu::particle_array_gpu<input_data_type> particle_array_type;

    /**
     * Creates a sampler for the given particle group and particle array.
     *
     * @param group     particle group providing the index list for the sample
     * @param array     particle array containing the data to be sampled
     * @param dim       CUDA configuration for launching kernels
     */
    phase_space_sampler_gpu_typed(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<mdsim::gpu::particle_array_gpu_base> array
      , cuda::config const& dim
    )
      : particle_group_(group), array_(mdsim::gpu::particle_array_gpu<input_data_type>::cast(array))
      , dim_(dim)
    {}

    /**
     * Static wrapper to directly construct a shared_ptr.
     */
    static std::shared_ptr<phase_space_sampler_gpu_typed> create(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<mdsim::gpu::particle_array_gpu_base> array
      , cuda::config const& dim
    )
    {
        return std::make_shared<phase_space_sampler_gpu_typed>(group, array, dim);
    }

    /**
     * acquire a sample
     *
     * copies the GPU data to a GPU sample using a CUDA kernel
     */
    virtual std::shared_ptr<sample_base> acquire()
    {
        if (!(cache_observer_ == array_->data()) || !(group_observer_ == particle_group_->ordered())) {
            auto const& data = read_cache(array_->data());
            auto const& group = read_cache(particle_group_->ordered());

            sample_ = std::make_shared<sample_type>(group.size());
            auto& sample_data = sample_->data();

            try {
                sample_wrapper::kernel.input.bind(data);
                cuda::configure(dim_.grid, dim_.block);
                sample_wrapper::kernel.sample(
                    &*group.begin()
                  , &*sample_data.begin()
                  , group.size()
                );
            }
            catch (cuda::error const&)
            {
                LOG_ERROR("failed to sample particle data on GPU");
                throw;
            }

            cache_observer_ = array_->data();
            group_observer_ = particle_group_->ordered();
        }
        return sample_;
    }

    /**
     * copies the data of a sample to the GPU particle array using a CUDA kernel
     */
    virtual void set(std::shared_ptr<sample_base const> sample_)
    {
        if(!sample_->gpu() || sample_->type() != typeid(sample_data_type)) {
            throw std::runtime_error("invalid sample data type");
        }
        auto sample = std::static_pointer_cast<sample_type const>(sample_);

        auto const& group = read_cache(particle_group_->ordered());

        auto data = make_cache_mutable(array_->mutable_data());
        try {
            sample_wrapper::kernel.input.bind(sample->data());
            cuda::configure(dim_.grid, dim_.block);
            sample_wrapper::kernel.set(
                &*group.begin()
              , data->data()
              , group.size()
            );
        }
        catch (cuda::error const&)
        {
            LOG_ERROR("failed to set particle data on GPU");
            throw;
        }
    }

    /**
     * returns a lua slot function to be used to acquire a gpu sample
     */
    virtual luaponte::object acquire_lua(lua_State* L, std::shared_ptr<phase_space_sampler_gpu> self)
    {
        std::function<std::shared_ptr<sample_type const>()> fn = [self]() -> std::shared_ptr<sample_type const>
        {
            return std::static_pointer_cast<sample_type const>(std::static_pointer_cast<phase_space_sampler_gpu_typed>(self)->acquire());
        };

        luaponte::default_converter<std::function<std::shared_ptr<sample_type const>()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }

    /**
     * returns a lua slot function to be used to directly acquire the data of a host sample
     */
    virtual luaponte::object data_lua(lua_State* L, std::shared_ptr<phase_space_sampler_gpu> self)
    {
        std::function<typename sample_type::array_type const&()> fn = [self]() -> typename sample_type::array_type const&
        {
            return std::static_pointer_cast<sample_type const>(std::static_pointer_cast<phase_space_sampler_gpu_typed>(self)->acquire())->data();
        };

        luaponte::default_converter<std::function<typename sample_type::array_type const&()>>().apply(L, fn);
        luaponte::object result(luaponte::from_stack(L, -1));
        lua_pop(L, 1);
        return result;
    }

    /**
     * wrapper to export the set member to lua
     */
    virtual void set_lua(luaponte::object sample)
    {
        set(luaponte::object_cast<std::shared_ptr<sample_type const>>(sample));
    }

protected:
    /** particle group containing the data indices to be sampled */
    std::shared_ptr<particle_group_type> particle_group_;
    /** particle array containing the data to be sampled */
    std::shared_ptr<particle_array_type> array_;
    /** cache observer for the array data */
    cache<> cache_observer_;
    /** cached GPU sample */
    std::shared_ptr<sample_type> sample_;
    /** cache observer for the index list */
    cache<> group_observer_;
    /** CUDA configuration for launching kernels */
    cuda::config const dim_;
};

/**
 * associative maps from typeid's to typed phase space sampler create functions for dimension 2 and 3
 *
 * used to create the correct phase space sampler based on the typeid of a particle array and the dimension of
 * the particle instance
 */
static const std::unordered_map<
    halmd::mdsim::gpu::ValueType
  , std::function<std::shared_ptr<phase_space_sampler_gpu>(
        std::shared_ptr<mdsim::gpu::particle_group>
      , std::shared_ptr<mdsim::gpu::particle_array_gpu_base>
      , cuda::config const&
    )>
> phase_space_sampler_gpu_typed_create_map[] = {
    { // dimension 2
        { halmd::mdsim::gpu::ValueType::FLOAT, phase_space_sampler_gpu_typed<1, float>::create }
      , { halmd::mdsim::gpu::ValueType::FLOAT2, phase_space_sampler_gpu_typed<2, float2>::create }
      , { halmd::mdsim::gpu::ValueType::FLOAT4, phase_space_sampler_gpu_typed<2, float4>::create }
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
      , { halmd::mdsim::gpu::ValueType::DSFLOAT, phase_space_sampler_gpu_typed<1, dsfloat, float>::create }
      , { halmd::mdsim::gpu::ValueType::DSFLOAT2, phase_space_sampler_gpu_typed<2, fixed_vector<dsfloat, 2>, float2>::create }
      , { halmd::mdsim::gpu::ValueType::DSFLOAT4, phase_space_sampler_gpu_typed<2, fixed_vector<dsfloat, 4>, float4>::create }
#endif
      , { halmd::mdsim::gpu::ValueType::INT, phase_space_sampler_gpu_typed<1, int>::create }
      , { halmd::mdsim::gpu::ValueType::INT2, phase_space_sampler_gpu_typed<2, int2>::create }
      , { halmd::mdsim::gpu::ValueType::INT4, phase_space_sampler_gpu_typed<2, int4>::create }

      , { halmd::mdsim::gpu::ValueType::UINT, phase_space_sampler_gpu_typed<1, unsigned int>::create }
      , { halmd::mdsim::gpu::ValueType::UINT2, phase_space_sampler_gpu_typed<2, uint2>::create }
      , { halmd::mdsim::gpu::ValueType::UINT4, phase_space_sampler_gpu_typed<2, uint4>::create }
    }
  , { // dimension 3
        { halmd::mdsim::gpu::ValueType::FLOAT, phase_space_sampler_gpu_typed<1, float>::create }
      , { halmd::mdsim::gpu::ValueType::FLOAT2, phase_space_sampler_gpu_typed<2, float2>::create }
      , { halmd::mdsim::gpu::ValueType::FLOAT4, phase_space_sampler_gpu_typed<3, float4>::create }

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
      , { halmd::mdsim::gpu::ValueType::DSFLOAT, phase_space_sampler_gpu_typed<1, dsfloat, float>::create }
      , { halmd::mdsim::gpu::ValueType::DSFLOAT2, phase_space_sampler_gpu_typed<2, fixed_vector<dsfloat, 2>, float2>::create }
      , { halmd::mdsim::gpu::ValueType::DSFLOAT4, phase_space_sampler_gpu_typed<3, fixed_vector<dsfloat, 4>, float4>::create }
#endif

      , { halmd::mdsim::gpu::ValueType::INT, phase_space_sampler_gpu_typed<1, int>::create }
      , { halmd::mdsim::gpu::ValueType::INT2, phase_space_sampler_gpu_typed<2, int2>::create }
      , { halmd::mdsim::gpu::ValueType::INT4, phase_space_sampler_gpu_typed<3, int4>::create }

      , { halmd::mdsim::gpu::ValueType::UINT, phase_space_sampler_gpu_typed<1, unsigned int>::create }
      , { halmd::mdsim::gpu::ValueType::UINT2, phase_space_sampler_gpu_typed<2, uint2>::create }
      , { halmd::mdsim::gpu::ValueType::UINT4, phase_space_sampler_gpu_typed<3, uint4>::create }
    }
};

/**
 * specialized phase_space sampler for gpu position data
 *
 * does the same as the generic GPU sampler, but additionally reduces/extends periodic positions
 */
template<int dimension, typename data_type>
class phase_space_sampler_gpu_position
  : public phase_space_sampler_gpu_typed<dimension, data_type, float4>
{
public:
    typedef samples::sample<dimension, float4> sample_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle_array_gpu<data_type> position_array_type;
    typedef mdsim::gpu::particle_array_gpu<typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type> image_array_type;

    /**
     * Creates a GPU sampler for position data.
     *
     * @param group             particle group providing the index list for the sample
     * @param box               box for periodic border conditions
     * @param position_array    particle array containing the position data
     * @param image_array       particle array containing the image data
     * @param dim       CUDA configuration for launching kernels
     */
    phase_space_sampler_gpu_position(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<mdsim::gpu::particle_array_gpu_base> position_array
      , std::shared_ptr<mdsim::gpu::particle_array_gpu_base> image_array
      , cuda::config const& dim
    )
      : phase_space_sampler_gpu_typed<dimension, data_type, float4>(group, position_array, dim)
      , box_(box)
      , image_array_(mdsim::gpu::particle_array_gpu<typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type>::cast(image_array))
    {}

    /**
     * Static wrapper to directly construct a shared_ptr.
     */
    static std::shared_ptr<phase_space_sampler_gpu_position> create(
        std::shared_ptr<particle_group_type> group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<mdsim::gpu::particle_array_gpu_base> position_array
      , std::shared_ptr<mdsim::gpu::particle_array_gpu_base> image_array
      , cuda::config const& dim
    )
    {
        return std::make_shared<phase_space_sampler_gpu_position>(group, box, position_array, image_array, dim);
    }

    /**
     * acquire a sample
     *
     * copies the GPU data to a GPU sample using a CUDA kernel
     */
    virtual std::shared_ptr<sample_base> acquire()
    {
        if (!(this->group_observer_ == this->particle_group_->ordered())
            || !(this->cache_observer_ == this->array_->data())
            || !(image_observer_ == image_array_->data())) {
            auto const& group = read_cache(this->particle_group_->ordered());
            auto const& particle_position = read_cache(this->array_->data());
            auto const& particle_image = read_cache(image_array_->data());
            this->group_observer_ = this->particle_group_->ordered();
            this->cache_observer_ = this->array_->data();
            image_observer_ = image_array_->data();

            this->sample_ = std::make_shared<sample_type>(group.size());

            try {
                phase_space_wrapper<dimension>::kernel.r.bind(particle_position);
                phase_space_wrapper<dimension>::kernel.image.bind(particle_image);

                cuda::configure(this->dim_.grid, this->dim_.block);
                phase_space_wrapper<dimension>::kernel.sample_position(
                    &*group.begin()
                  , this->sample_->data()
                  , static_cast<fixed_vector<float, dimension>>(box_->length())
                  , group.size()
                );
            }
            catch (cuda::error const&)
            {
                LOG_ERROR("failed to sample particle positions on GPU");
                throw;
            }

        }
        return this->sample_;
    }

    /**
     * copies the data of a sample to the GPU particle array using a CUDA kernel
     */
    virtual void set(std::shared_ptr<sample_base const> sample)
    {
        // set position data the same as any other kind of data
        phase_space_sampler_gpu_typed<dimension, data_type, float4>::set(sample);

        // reduce positions on GPU
        auto position = make_cache_mutable(this->array_->mutable_data());
        auto image = make_cache_mutable(this->image_array_->mutable_data());
        auto const& group = read_cache(this->particle_group_->ordered());
        try {
            phase_space_wrapper<dimension>::kernel.r.bind(*position);
            cuda::configure(this->dim_.grid, this->dim_.block);
            phase_space_wrapper<dimension>::kernel.reduce_periodic(
                &*group.begin()
              , position->data() // TODO: is this correct for dsfloats?
              , &*image->begin()
              , static_cast<fixed_vector<float, dimension>>(box_->length())
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
    /** box for periodic boundary conditions */
    std::shared_ptr<box_type const> box_;
    /** particle array containing the image data */
    std::shared_ptr<image_array_type> image_array_;
    /** cache observer for the image data */
    cache<> image_observer_;
};


template <typename phase_space_type>
static int wrap_dimension(phase_space_type const&)
{
    return phase_space_type::particle_type::vector_type::static_size;
}

template <int dimension, typename float_type>
phase_space<dimension, float_type>::phase_space(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<particle_group_type> particle_group
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , particle_group_(particle_group)
  , box_(box)
  , logger_(logger)
{}

/**
 * Get phase space sampler implementation.
 */
template <int dimension, typename float_type>
std::shared_ptr<phase_space_sampler_gpu>
phase_space<dimension, float_type>::get_sampler_gpu(std::string const& name)
{
    auto it = gpu_samplers_.find(name);

    if (it != gpu_samplers_.end()) {
        return it->second;
    } else {
        auto array = particle_->get_gpu_array(name);
        if(!name.compare("position")) {
            return (gpu_samplers_[name] = phase_space_sampler_gpu_position<dimension, typename mdsim::gpu::particle<dimension, float_type>::gpu_hp_vector_type>::create
                    (particle_group_, box_, array, particle_->get_gpu_array("image"), particle_->dim()));
        } else {
            auto it = phase_space_sampler_gpu_typed_create_map[dimension-2].find(array->value_type());
            if(it == phase_space_sampler_gpu_typed_create_map[dimension-2].end()) {
                throw std::runtime_error("invalid sample type");
            }
            return (gpu_samplers_[name] = it->second(particle_group_, array, particle_->dim()));
        }
    }
}

template <int dimension, typename float_type>
std::shared_ptr<phase_space_sampler_host>
phase_space<dimension, float_type>::get_sampler_host(std::string const& name)
{
    auto it = host_samplers_.find(name);

    if (it != host_samplers_.end()) {
        return it->second;
    } else {
        auto array = particle_->get_host_array(name);
        if(!name.compare("position")) {
            return (host_samplers_[name] = phase_space_sampler_position<dimension, float_type, float>::create
                    (particle_group_, box_, array, particle_->get_host_array("image"), particle_->dim(), host_cache_));
        } else {
            auto it = phase_space_sampler_typed_create_map.find(array->type());
            if(it == phase_space_sampler_typed_create_map.end()) {
                throw std::runtime_error("invalid sample type");
            }
            return (host_samplers_[name] = it->second(particle_group_, array, particle_->dim().threads(), host_cache_));
        }
    }
}

template <typename phase_space_type>
static luaponte::object wrap_acquire(lua_State* L, std::shared_ptr<phase_space_type> self, std::string const& name, bool gpu)
{
    if (gpu) {
        auto sampler = self->get_sampler_gpu(name);
        return sampler->acquire_lua(L, sampler);
    } else {
        auto sampler = self->get_sampler_host(name);
        return sampler->acquire_lua(L, sampler);
    }
}

template <typename phase_space_type>
static void wrap_set(std::shared_ptr<phase_space_type> self, std::string const& name, luaponte::object sample, bool gpu)
{
    if (gpu) {
        self->get_sampler_gpu(name)->set_lua(sample);
    } else {
        self->get_sampler_host(name)->set_lua(sample);
    }
}

template <typename phase_space_type>
static luaponte::object wrap_data(lua_State* L, std::shared_ptr<phase_space_type> self, std::string const& name)
{
    auto sampler = self->get_sampler_host(name);
    return sampler->data_lua(L, sampler);
}

template <typename phase_space_type>
static luaponte::object wrap_gpu_data(lua_State* L, std::shared_ptr<phase_space_type> self, std::string const& name)
{
    auto sampler = self->get_sampler_gpu(name);
    return sampler->data_lua(L, sampler);
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
                .def("gpu_data", &wrap_gpu_data<phase_space>)
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

          , def("phase_space", &std::make_shared<phase_space
                , std::shared_ptr<particle_type>
                , std::shared_ptr<particle_group_type>
                , std::shared_ptr<box_type const>
                , std::shared_ptr<logger>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_phase_space(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    phase_space<3, dsfloat>::luaopen(L);
    phase_space<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class phase_space<3, float>;
template class phase_space<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class phase_space<3, dsfloat>;
template class phase_space<2, dsfloat>;
#endif

} // namespace gpu
} // namespace observables
} // namespace halmd
