/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_MORSE_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_MORSE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/morse_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

/**
 * define Morse potential and parameters
 */
template <typename float_type>
class morse
{
public:
    typedef morse_kernel::morse gpu_potential_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef logger logger_type;

    morse(
        unsigned int ntype1
      , unsigned int ntype2
      , matrix_type const& cutoff
      , matrix_type const& epsilon
      , matrix_type const& sigma
      , matrix_type const& r_min
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /** bind textures before kernel invocation */
    void bind_textures() const
    {
        morse_wrapper::param.bind(g_param_);
        morse_wrapper::rr_cut.bind(g_rr_cut_);
    }

    matrix_type const& r_cut() const
    {
        return r_cut_;
    }

    float_type r_cut(unsigned a, unsigned b) const
    {
        return r_cut_(a, b);
    }

    float_type rr_cut(unsigned a, unsigned b) const
    {
        return rr_cut_(a, b);
    }

    matrix_type const& r_cut_sigma() const
    {
        return r_cut_sigma_;
    }

    matrix_type const& epsilon() const
    {
        return epsilon_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

    matrix_type const& r_min_sigma() const
    {
        return r_min_sigma_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** depths of potential well in MD units */
    matrix_type epsilon_;
    /** width of potential well in MD units */
    matrix_type sigma_;
    /** position of potential well in units of sigma */
    matrix_type r_min_sigma_;
    /** cutoff radius in units of sigma */
    matrix_type r_cut_sigma_;
    /** cutoff radius in MD units */
    matrix_type r_cut_;
    /** square of cutoff radius */
    matrix_type rr_cut_;
    /** potential energy at cutoff length in MD units */
    matrix_type en_cut_;
    /** potential parameters at CUDA device */
    cuda::vector<float4> g_param_;
    /** squared cutoff radius at CUDA device */
    cuda::vector<float> g_rr_cut_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;
};

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_MORSE_HPP */
