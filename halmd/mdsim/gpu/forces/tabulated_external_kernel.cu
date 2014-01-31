/*
 * Copyright © 2013 Nicolas Höft
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

#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/forces/interpolation/cubic_hermite.hpp>
#include <halmd/mdsim/forces/interpolation/linear.hpp>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/forces/tabulated_external_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {
namespace tabulated_external_kernel {

/**
 * Compute pair forces, potential energy, and stress tensor for all particles
 */
template <
    bool do_aux               //< compute auxiliary variables in addition to force
  , typename vector_type
  , typename gpu_vector_type
  , typename force_interpolation_type
  , typename virial_interpolation_type
>
__global__ void compute(
    float4 const* g_r
  , gpu_vector_type* g_f
  , float* g_en_pot
  , float* g_stress_pot
  , unsigned int ntype
  , vector_type box_length
  , float const* g_force_coefficients
  , float const* g_virial_coefficients
  , force_interpolation_type force_interpolation
  , virial_interpolation_type virial_interpolation
  , bool force_zero
)
{
    enum { dimension = vector_type::static_size };
    typedef typename vector_type::value_type value_type;
    typedef typename force_interpolation_type::index_type grid_index_type;
    typedef typename type_traits<dimension, float>::stress_tensor_type stress_tensor_type;
/*#ifdef USE_FORCE_DSFUN
    typedef fixed_vector<dsfloat, dimension> force_vector_type;
#else*/
    typedef vector_type force_vector_type;
//#endif

    unsigned int const i = GTID;

    // load particle associated with this thread
    unsigned int type;
    vector_type r;
    tie(r, type) <<= g_r[i];

    // contribution to potential energy
    float en_pot_ = 0;
    // contribution to stress tensor
    stress_tensor_type stress_pot = 0;

    // force sum
    force_vector_type f = 0;

    // enforce periodic boundary conditions
    box_kernel::reduce_periodic(r, box_length);

    // apply interpolation
    // http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/index.jsp?topic=/com.ibm.xlcpp8a.doc/language/ref/keyword_template_qualifier.htm
    tie(en_pot_, f) = force_interpolation.template operator()<force_vector_type>(r, g_force_coefficients);

    if (do_aux) {
        // contribution to stress tensor from this particle
        float virial;
        force_vector_type dummy;
        tie(virial, dummy) = virial_interpolation.template operator()<force_vector_type>(r, g_virial_coefficients);
        for (int d = 0; d < dimension; ++d) {
            stress_pot[d] = virial/dimension;
        }
    }

    // add old force and auxiliary variables if not zero
    if (!force_zero) {
        f += static_cast<vector_type>(g_f[i]);
        if(do_aux) {
            en_pot_ += g_en_pot[i];
            stress_pot += read_stress_tensor<stress_tensor_type>(g_stress_pot + i, GTDIM);
        }
    }

    // write results to global memory
    g_f[i] = static_cast<vector_type>(f);
    if (do_aux) {
        g_en_pot[i] = en_pot_;
        write_stress_tensor(g_stress_pot + i, stress_pot, GTDIM);
    }
}

} // namespace tabulated_external_kernel

template <int dimension, typename interpolation_type, typename virial_interpolation_type>
tabulated_external_wrapper<dimension, interpolation_type, virial_interpolation_type> const
tabulated_external_wrapper<dimension, interpolation_type, virial_interpolation_type>::kernel = {
    tabulated_external_kernel::compute<false, fixed_vector<float, dimension> >
  , tabulated_external_kernel::compute<true, fixed_vector<float, dimension> >
};

using namespace halmd::mdsim::forces::interpolation;
// explicit instantiation
template class tabulated_external_wrapper<3, cubic_hermite<3, float>, linear<3, float> >;
template class tabulated_external_wrapper<2, cubic_hermite<2, float>, linear<2, float> >;
template class tabulated_external_wrapper<3, cubic_hermite<3, double>, linear<3, float> >;
template class tabulated_external_wrapper<2, cubic_hermite<2, double>, linear<2, float> >;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd

