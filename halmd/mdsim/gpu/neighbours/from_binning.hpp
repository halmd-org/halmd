/*
 * Copyright © 2008-2015  Felix Höfling
 * Copyright © 2015       Nicolas Höft
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_HPP
#define HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/binning.hpp>
#include <halmd/mdsim/gpu/max_displacement.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/profiler.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <memory>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

template <int dimension, typename float_type>
class from_binning
  : public gpu::neighbour
{
private:
    typedef gpu::neighbour _Base;

public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef mdsim::box<dimension> box_type;
    typedef gpu::binning<dimension, float_type> binning_type;
    typedef max_displacement<dimension, float_type> displacement_type;
    struct defaults;

    enum algorithm
    {
        shared_mem = 1
      , naive      = 2
    };

    typedef _Base::array_type array_type;

    static void luaopen(lua_State* L);

    from_binning(
        std::pair<std::shared_ptr<particle_type const>, std::shared_ptr<particle_type const>> particle
      , std::pair<std::shared_ptr<binning_type>, std::shared_ptr<binning_type>> binning
      , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>> displacement
      , std::shared_ptr<box_type const> box
      , matrix_type const& r_cut
      , double skin
      , double cell_occupancy = defaults::occupancy()
      , algorithm preferred_algorithm = shared_mem
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    connection on_prepend_update(std::function<void ()> const& slot)
    {
        return on_prepend_update_.connect(slot);
    }

    connection on_append_update(std::function<void ()> const& slot)
    {
        return on_append_update_.connect(slot);
    }

    //! returns neighbour list skin in MD units
    float_type r_skin() const
    {
        return r_skin_;
    }

    //! returns neighbour list skin in MD units
    float_type cell_occupancy() const
    {
        return nu_cell_;
    }

    /**
     * neighbour lists
     */
    virtual cache<array_type> const& g_neighbour();

    /**
     * number of placeholders per neighbour list
     */
    virtual unsigned int size() const
    {
        return size_;
    }

    /**
     * neighbour list stride
     */
    virtual unsigned int stride() const
    {
        return stride_;
    }

    //! returns true if the binning modules are compatible with the neighbour list module
    static bool is_binning_compatible(
        std::shared_ptr<binning_type const> binning1
      , std::shared_ptr<binning_type const> binning2
    );

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    typedef typename binning_type::array_type cell_array_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type update;
    };

    void update();
    void set_occupancy(double occupancy);

    std::shared_ptr<particle_type const> particle1_;
    std::shared_ptr<particle_type const> particle2_;
    std::shared_ptr<binning_type> binning1_;
    std::shared_ptr<binning_type> binning2_;
    std::shared_ptr<displacement_type> displacement1_;
    std::shared_ptr<displacement_type> displacement2_;
    std::shared_ptr<box_type const> box_;
    std::shared_ptr<logger> logger_;

    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** maximum cutoff length */
    float_type r_cut_max_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    cuda::vector<float_type> g_rr_cut_skin_;
    /** FIXME average desired cell occupancy */
    float_type nu_cell_;
    /** preferred algorithm for update */
    algorithm preferred_algorithm_;
    /** neighbour lists */
    cache<array_type> g_neighbour_;
    /** cache observer for neighbour list update */
    std::tuple<cache<>, cache<>> neighbour_cache_;
    /** number of placeholders per neighbour list */
    unsigned int size_;
    /** neighbour list stride */
    unsigned int stride_;

    /** CUDA device properties */
    cuda::device::properties device_properties_;

    /** profiling runtime accumulators */
    runtime runtime_;
    /** signal emitted before neighbour list update */
    signal<void ()> on_prepend_update_;
    /** signal emitted after neighbour list update */
    signal<void ()> on_append_update_;
};

template <int dimension, typename float_type>
struct from_binning<dimension, float_type>::defaults
{
    static float_type occupancy();
};

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_HPP */
