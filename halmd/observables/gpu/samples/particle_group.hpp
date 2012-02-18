/*
 * Copyright © 2012  Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_OBSERVABLES_GPU_SAMPLES_PARTICLE_GROUP_HPP
#define HALMD_OBSERVABLES_GPU_SAMPLES_PARTICLE_GROUP_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <utility>

#include <halmd/mdsim/gpu/particle.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

/**
 * A particle group represents a subset of particles, which is defined here by
 * an instance of gpu::particle together with either a range of tags or by
 * selecting all.
 *
 * A tag range is a contiguous range of particle tags, specified in terms of
 * begin and end tags in analogy to iterator ranges, the particle with tag
 * 'begin' is included, while tag 'end' is not.
 *
 * The group represents a fixed order of the particles according to their tags
 * and starts with the smallest tag in the set.
 *
 */

template <int dimension, typename float_type>
class particle_group
{
public:
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef cuda::vector<unsigned int>::const_pointer gpu_map_iterator; // FIXME vector<T>::const_iterator

    static void luaopen(lua_State* L);

    particle_group() {}

    //! returns underlying particle instance
    virtual boost::shared_ptr<particle_type const> particle() const = 0;

    /**
     * returns iterator to an index array in GPU memory mapping particle tags
     * to array indices in gpu::particle
     */
    virtual gpu_map_iterator g_map() const = 0;
    /**
     * returns iterator to an index array in host memory mapping particle tags
     * to array indices in gpu::particle
     *
     * Method can not be const since it copies to internal buffer in host
     * memory.
     */
    virtual unsigned int const* h_map() = 0;

    //! returns the size of the group, i.e., the number of particles
    virtual unsigned int size() const = 0;

    //! returns true if the group is the empty set
    bool empty() const
    {
        return size() == 0;
    }
};

template <int dimension, typename float_type>
class particle_group_all
  : public particle_group<dimension, float_type>
{
public:
    typedef particle_group<dimension, float_type> _Base;
    typedef typename _Base::particle_type particle_type;
    typedef typename _Base::gpu_map_iterator gpu_map_iterator;

    static void luaopen(lua_State* L);

    particle_group_all(
        boost::shared_ptr<particle_type const> particle
    );

    virtual boost::shared_ptr<particle_type const> particle() const
    {
        return particle_;
    }

    virtual gpu_map_iterator g_map() const
    {
        return particle_->g_reverse_tag.data();
    }

    virtual unsigned int const* h_map();

    //! returns size of the group, i.e., the number of particles
    virtual unsigned int size() const
    {
        return particle_->nbox;
    }

private:
    /** gpu::particle instance */
    boost::shared_ptr<particle_type const> particle_;
    /** pre-allocate page-locked memory for reverse particle tags */
    cuda::host::vector<unsigned int> h_reverse_tag_;
};

template <int dimension, typename float_type>
class particle_group_from_range
  : public particle_group<dimension, float_type>
{
public:
    typedef particle_group<dimension, float_type> _Base;
    typedef typename _Base::particle_type particle_type;
    typedef typename _Base::gpu_map_iterator gpu_map_iterator;

    static void luaopen(lua_State* L);

    particle_group_from_range(
        boost::shared_ptr<particle_type const> particle
      , unsigned int begin
      , unsigned int end
    );

    virtual boost::shared_ptr<particle_type const> particle() const
    {
        return particle_;
    }

    virtual gpu_map_iterator g_map() const
    {
        return particle_->g_reverse_tag.data() + begin_;
    }

    virtual unsigned int const* h_map();

    //! returns size of the group, i.e., the number of particles
    virtual unsigned int size() const
    {
        return end_ - begin_;
    }

private:
    /** gpu::particle instance */
    boost::shared_ptr<particle_type const> particle_;
    /** tag range [begin, end) */
    unsigned int begin_;
    unsigned int end_;

    /** pre-allocate page-locked memory for reverse particle tags */
    cuda::host::vector<unsigned int> h_reverse_tag_;
};
} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_SAMPLES_PARTICLE_GROUP_HPP */
