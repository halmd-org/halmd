/*
 * Copyright © 2008-2010  Peter Colberg
 * Copyright © 2015       Nicolas Höft
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

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/gpu/phase_space_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include "phase_space_kernel.hpp"

using namespace halmd::mdsim::gpu; //< namespace box_kernel

namespace halmd {
namespace observables {
namespace gpu {
namespace phase_space_kernel {

/** positions, types */
texture<float4> r_;

/** minimum image vectors */
template<int dimension>
struct image
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename phase_space_wrapper<dimension>::coalesced_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> image<dimension>::type image<dimension>::tex_;

/**
 * sample phase space for all particle of a single species
 */
template <typename vector_type, typename T>
__global__ void sample_position(
    unsigned int const* g_reverse_id
  , T* g_r
  , vector_type box_length
  , unsigned int npart
)
{
    enum { dimension = vector_type::static_size };
    typedef typename phase_space_wrapper<dimension>::coalesced_vector_type coalesced_vector_type;

    if (GTID < npart) {
        // permutation index
        uint const rid = g_reverse_id[GTID];
        // fetch particle from texture caches
        unsigned int type;
        vector_type r;
        tie(r, type) <<= tex1Dfetch(r_, rid);
        // extend particle positions in periodic box
        vector_type img = tex1Dfetch(image<dimension>::tex_, rid);
        box_kernel::extend_periodic(r, img, box_length);
        // store particle in global memory
        g_r[GTID] <<= tie(r, type);
    }
}

/**
 * shift particle positions to range (-L/2, L/2)
 */
template <typename vector_type, typename coalesced_vector_type>
__global__ void reduce_periodic(
    unsigned int const* g_reverse_id
  , float4* g_r
  , coalesced_vector_type* g_image
  , vector_type box_length
  , unsigned int npart
)
{
    enum { dimension = vector_type::static_size };

    if (GTID < npart) {
        unsigned int rid = g_reverse_id[GTID];
        vector_type r;
        unsigned int type;
        tie(r, type) <<= tex1Dfetch(r_, rid);

        vector_type image = box_kernel::reduce_periodic(r, box_length);

        g_image[rid] = image;
        g_r[rid] <<= tie(r, type);
    }
}

} // namespace phase_space_kernel

template <int dimension>
phase_space_wrapper<dimension> const phase_space_wrapper<dimension>::kernel = {
    phase_space_kernel::r_
  , phase_space_kernel::image<dimension>::tex_
  , phase_space_kernel::sample_position<fixed_vector<float, dimension> >
  , phase_space_kernel::reduce_periodic<fixed_vector<float, dimension> >
};

template class phase_space_wrapper<3>;
template class phase_space_wrapper<2>;

namespace phase_space_sample_kernel {

/** minimum image vectors */
template<typename T>
struct input
{
    // instantiate a separate texture for each vector type
    typedef texture<T> type;
    static type tex_;
};
// instantiate static members
template<typename T> input<T>::type input<T>::tex_;

template<typename T, typename U>
struct converter
{
    static __device__ T const& convert (U const& u)
    {
        return u;
    }
};

template<size_t dimension, typename U>
struct converter<fixed_vector<dsfloat, dimension>, U>
{
    static __device__ tuple<U,U> convert (U const& u)
    {
        return make_tuple(u, U());
    }
};

template<typename U>
struct converter<dsfloat, U> : converter<fixed_vector<dsfloat, 1>, U> {};

template <typename ptr_type, typename vector_type, typename T>
__global__ void sample(
        unsigned int const* g_reverse_id
        , ptr_type data
        , unsigned int npart
) {
    if (GTID < npart) {
        // permutation index
        uint const rid = g_reverse_id[GTID];
        // fetch particle data from texture caches
        data[GTID] = converter<vector_type, T>::convert(tex1Dfetch(input<T>::tex_, rid));
    }
}

} // namespace phase_space_sample_kernel

template <typename input_data_type, typename sample_data_type>
phase_space_sample_wrapper<input_data_type, sample_data_type> const phase_space_sample_wrapper<input_data_type, sample_data_type>::kernel = {
        phase_space_sample_kernel::input<sample_data_type>::tex_,
        phase_space_sample_kernel::sample<sample_data_type*, sample_data_type, sample_data_type>,
        phase_space_sample_kernel::sample<ptr_type, input_data_type, sample_data_type>
};

template class phase_space_sample_wrapper<float>;
template class phase_space_sample_wrapper<float2>;
template class phase_space_sample_wrapper<float4>;

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class phase_space_sample_wrapper<dsfloat, float>;
template class phase_space_sample_wrapper<fixed_vector<dsfloat, 2>, float2>;
template class phase_space_sample_wrapper<fixed_vector<dsfloat, 4>, float4>;
#endif

template class phase_space_sample_wrapper<int>;
template class phase_space_sample_wrapper<int2>;
template class phase_space_sample_wrapper<int4>;

template class phase_space_sample_wrapper<unsigned int>;
template class phase_space_sample_wrapper<uint2>;
template class phase_space_sample_wrapper<uint4>;

} // namespace gpu
} // namespace observables
} // namespace halmd
