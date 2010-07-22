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

#ifndef HALMD_MDSIM_GPU_NEIGHBOUR_HPP
#define HALMD_MDSIM_GPU_NEIGHBOUR_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/gpu/box.hpp>
#include <halmd/mdsim/gpu/neighbour_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/neighbour.hpp>
// #include <halmd/numeric/gpu/blas/vector.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

// forward declaration
template <int dimension, typename float_type>
class force;

namespace sort
{
template <int dimension, typename float_type>
class hilbert;
}

template <int dimension, typename float_type>
class neighbour
  : public mdsim::neighbour<dimension>
{
public:
    // module definitions
    typedef neighbour _Self;
    typedef mdsim::neighbour<dimension> _Base;
    static void options(po::options_description& desc);
    static void depends();
    static void select(po::options const& vm) {}

    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef gpu::force<dimension, float_type> force_type;
    typedef typename force_type::matrix_type matrix_type;
    typedef gpu::box<dimension> box_type;

    typedef typename neighbour_wrapper<dimension>::displacement_impl_type displacement_impl_type;

//     typedef typename particle_type::neighbour_list cell_list;
//     typedef boost::multi_array<cell_list, dimension> cell_lists;
    typedef fixed_vector<unsigned int, dimension> cell_size_type;
    typedef fixed_vector<int, dimension> cell_diff_type;

    shared_ptr<particle_type> particle;
    shared_ptr<force_type> force;
    shared_ptr<box_type> box;

    cuda::config dim_reduce;
    displacement_impl_type const displacement_impl;

    neighbour(modules::factory& factory, po::options const& vm);
    virtual ~neighbour() {}
    void update();
    bool check();

protected:
    friend class sort::hilbert<dimension, float_type>;

    static displacement_impl_type get_displacement_impl(int threads);
    void update_cells();
    void update_neighbours();
//     void update_cell_neighbours(cell_size_type const& i);
//     template <bool same_cell>
//     void compute_cell_neighbours(size_t i, cell_list& c);

    /** neighbour list skin in MD units */
    float_type r_skin_;
    /** half neighbour list skin */
    float_type rr_skin_half_;
    /** (cutoff lengths + neighbour list skin)² */
    matrix_type rr_cut_skin_;
    /** (cutoff lengths + neighbour list skin)² */
    cuda::vector<float_type> g_rr_cut_skin_;
    /** average desired cell occupancy */
    float_type nu_cell_;
    /** number of cells per dimension */
    cell_size_type ncell_;
    /** number of placeholders per cell */
    size_t cell_size_;
    /** cell edge lengths */
    vector_type cell_length_;
    /** CUDA cell kernel execution configuration */
    cuda::config dim_cell_;
    /** particle positions at last neighbour list update */
    cuda::vector<float4> g_r0_;
    /** block-reduced squared particle distances */
    cuda::vector<float> g_rr_;
    /** block-reduced squared particle distances */
    cuda::host::vector<float> h_rr_;
    /** cell lists in global device memory */
    cuda::vector<unsigned int> g_cell_;

    /** GPU radix sort */
    algorithm::gpu::radix_sort<unsigned int> sort_;
    /** cell indices for particles */
    cuda::vector<unsigned int> g_cell_index_;
    /** particle permutation */
    cuda::vector<unsigned int> g_cell_permutation_;
    /** cell offsets in sorted particle list */
    cuda::vector<unsigned int> g_cell_offset_;
};

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOUR_HPP */
