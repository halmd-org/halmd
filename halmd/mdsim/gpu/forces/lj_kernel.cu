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

#include <halmd/mdsim/gpu/forces/lj_kernel.hpp>
#include <halmd/mdsim/gpu/forces/pair_short_ranged_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/variant.cuh>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{
namespace lj_kernel
{

/** array of Lennard-Jones potential parameters for all combinations of particle types */
static texture<float4> param_;

/** define Lennard-Jones potential */
struct lj_potential
{
    template <typename float_type, typename param_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator() (float_type rr, param_type const& param) const
    {
        float_type rri = param[SIGMA2] / rr;
        float_type ri6 = rri * rri * rri;
        float_type fval = 48 * param[EPSILON] * rri * ri6 * (ri6 - 0.5f) / param[SIGMA2];
        float_type en_pot = 4 * param[EPSILON] * ri6 * (ri6 - 1) - param[EN_CUT];

        return make_tuple(fval, en_pot);
    }

    HALMD_GPU_ENABLED texture<float4> const& param() const
    {
        return param_;
    }
};

} // namespace lj_kernel

template <int dimension>
lj_wrapper<dimension> const lj_wrapper<dimension>::kernel = {
    lj_kernel::param_
};

template class lj_wrapper<3>;
template class lj_wrapper<2>;

template class pair_short_ranged_wrapper<3, lj_kernel::lj_potential>;
template class pair_short_ranged_wrapper<2, lj_kernel::lj_potential>;

}}} // namespace mdsim::gpu::forces

} // namespace halmd
