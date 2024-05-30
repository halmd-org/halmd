/*
 * Copyright © 2014-2015 Sutapa Roy
 * Copyright © 2014-2015 Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_PLANAR_WALL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_PLANAR_WALL_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/external/planar_wall_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {

/**
 * define planar_wall potential and parameters
 */
template <int dimension, typename float_type>
class planar_wall
{
public:
    typedef planar_wall_kernel::planar_wall<dimension> gpu_potential_type;

    typedef fixed_vector<float_type, dimension> vector_type;
    typedef boost::numeric::ublas::vector<float_type> scalar_container_type;
    typedef boost::numeric::ublas::vector<vector_type> vector_container_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_container_type;

    planar_wall(
        scalar_container_type const& offset
      , vector_container_type const& surface_normal
      , matrix_container_type const& epsilon
      , matrix_container_type const& sigma
      , matrix_container_type const& wetting
      , matrix_container_type const& cutoff
      , float_type smoothing
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /** bind textures before kernel invocation */
    void bind_textures() const
    {
        planar_wall_wrapper::param_geometry.bind(g_param_geometry_);
        planar_wall_wrapper::param_potential.bind(g_param_potential_);
    }

    scalar_container_type const& offset() const
    {
        return offset_;
    }

    vector_container_type const& surface_normal() const
    {
        return surface_normal_;
    }

    matrix_container_type const& epsilon() const
    {
        return epsilon_;
    }

    matrix_container_type const& sigma() const
    {
        return sigma_;
    }

    matrix_container_type const& wetting() const
    {
        return wetting_;
    }

    matrix_container_type const& cutoff() const
    {
        return cutoff_;
    }

    float_type smoothing()
    {
        return smoothing_;
    }

    unsigned int size() const
    {
        return offset_.size();
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** wall position in MD units */
    scalar_container_type offset_;
    /** wall normal vector in MD units */
    vector_container_type surface_normal_;
    /** interaction strength for wall potential in MD units */
    matrix_container_type epsilon_;
    /** interaction range for wall potential in MD units */
    matrix_container_type sigma_;
    /** wetting parameter for wall potential in MD units */
    matrix_container_type wetting_;
    /** cutoff length for wall potential in MD units */
    matrix_container_type cutoff_;
    /** smoothing parameters for wall potential in MD units */
    float_type smoothing_;

    /** potential parameters at CUDA device */
    cuda::vector<float4> g_param_geometry_;
    cuda::vector<float4> g_param_potential_;

    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace external
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_PLANAR_WALL_HPP */
