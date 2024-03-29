/*
 * Copyright © 2008-2023 Felix Höfling
 * Copyright © 2008-2010 Peter Colberg
 * Copyright © 2020      Jaslo Ziska
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/potentials/pair/mie_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * define Mie potential and parameters
 *
 * @f[ U(r) = C(m, n) \epsilon \left[ (r/\sigma)^{-m} - (r/\sigma)^{-n}) \right] @f]
 *
 * with prefactor @f$ C(m, n) = \frac{m}{m - n} \left(\frac{m}{n}\right)^{n/(m-n)}$f@.
 * The exponents @f$ m, n @f$ must be even and @f$ m > n @f$.
 */
template <typename float_type_>
class mie
{
public:
    typedef float_type_ float_type;
    typedef mie_kernel::mie gpu_potential_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef boost::numeric::ublas::matrix<unsigned> uint_matrix_type;

    mie(
        matrix_type const& epsilon
      , matrix_type const& sigma
      , uint_matrix_type const& index_m
      , uint_matrix_type const& index_n
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /** return gpu potential with texture */
    gpu_potential_type get_gpu_potential()
    {
        // FIXME: tex1Dfetch reads zero when texture is not recreated once in a while
        t_param_ = cuda::texture<float4>(g_param_);
        return gpu_potential_type(t_param_);
    }

    matrix_type const& epsilon() const
    {
        return epsilon_;
    }

    matrix_type const& sigma() const
    {
        return sigma_;
    }

    uint_matrix_type const& index_m() const
    {
        return index_m_;
    }

    uint_matrix_type const& index_n() const
    {
        return index_n_;
    }

    unsigned int size1() const
    {
        return epsilon_.size1();
    }

    unsigned int size2() const
    {
        return epsilon_.size2();
    }

    std::tuple<float_type, float_type> operator()(float_type rr, unsigned a, unsigned b) const
    {
        return mie_kernel::compute(rr
          , sigma_(a, b) * sigma_(a, b), epsilon_C_(a, b)
          , index_m_(a, b) / 2, index_n_(a, b) / 2
        );
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** potential well depths in MD units */
    matrix_type epsilon_;
    /** potential well depths times prefactor C(m,n) */
    matrix_type epsilon_C_;
    /** pair separation in MD units */
    matrix_type sigma_;
    /** power law index of repulsion */
    uint_matrix_type index_m_;
    /** power law index of attraction */
    uint_matrix_type index_n_;
    /** square of pair separation */
    matrix_type sigma2_;
    /** potential parameters at CUDA device */
    cuda::memory::device::vector<float4> g_param_;
    /** array of Lennard-Jones potential parameters for all combinations of particle types */
    cuda::texture<float4> t_param_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_HPP */
