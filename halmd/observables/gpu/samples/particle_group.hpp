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

#include <lua.hpp>
#include <utility>

#include <halmd/mdsim/gpu/particle.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

class range_tag : public std::pair<unsigned int, unsigned int> {};
struct all_tag {};
// struct unsorted_all_tag {};

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

    static void luaopen(lua_State* L);

    particle_group(
        boost::shared_ptr<particle_type /* FIXME const */> particle
      , range_tag range
    );

    particle_group(
        boost::shared_ptr<particle_type /* FIXME const */> particle
      , all_tag
    );

    boost::shared_ptr<particle_type /* FIXME const*/> particle() const
    {
        return particle_;
    }

    /**
     * returns a permutation array in GPU memory to map particle tags to array
     * indices in gpu::particle
     */
    unsigned int const* g_map() const
    {
        return particle_->g_reverse_tag.data() + begin_;
    }
    /**
     * returns a permutation array in host memory to map particle tags to array
     * indices in gpu::particle
     */
    unsigned int const* h_map() const;

    //! returns size of the group, i.e., the number of particles
    unsigned int size() const
    {
        return end_ - begin_;
    }

    //! returns true if group is the empty set
    bool empty() const
    {
        return begin_ == end_;
    }

private:
    /** gpu::particle instance */
    boost::shared_ptr<particle_type /* FIXME const */> particle_;
    /** tag range [begin, end) */
    unsigned int begin_;
    unsigned int end_;
};

} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_SAMPLES_PARTICLE_GROUP_HPP */
