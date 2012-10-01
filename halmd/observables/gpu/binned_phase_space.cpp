/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <functional>

#include <halmd/utility/gpu/device.hpp>
#include <halmd/observables/gpu/binned_phase_space.hpp>
#include <halmd/observables/gpu/binned_phase_space_kernel.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
binned_phase_space<dimension, float_type>::binned_phase_space(
    shared_ptr<sample_type> sample
  , shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , unsigned int threads
  , shared_ptr<logger_type> logger
)
  // module dependencies
  : sample_(sample)
  , clock_(clock)
  , logger_(logger)
  // initialise attributes
  , npart_(sample_->position_data().size()) // number of particles
  , dim_(device::validate(
        cuda::config((npart_ + threads - 1) / threads, threads)  // default CUDA configuration
    ))
  , sort_(npart_, threads)
  , g_cell_offset_(sample_->g_cell.size())
{
    LOG("using " << dim_.threads_per_block() << " CUDA threads in " << dim_.blocks_per_grid() << " blocks");

    vector_type box_length = static_cast<vector_type>(box->length());
    sample_->cell_length = element_div(box_length, static_cast<vector_type>(sample_->nbin));
    sample_->cell_origin = static_cast<vector_type>(box->origin());

    // allocate memory
    try {
        g_cell_index_.reserve(dim_.threads());
        g_cell_index_.resize(npart_);
        g_cell_permutation_.reserve(dim_.threads());
        g_cell_permutation_.resize(npart_);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate global device memory for phase space binning");
        throw;
    }

    // set up grid of cell positions axis-wise
    for (unsigned axis = 0; axis < dimension; ++axis) {
        double spacing = sample_->cell_length[axis];
        double offset = sample_->cell_origin[axis];
        position_[axis].reserve(sample_->nbin[axis]);
        for (unsigned i = 0; i < sample_->nbin[axis]; ++i) {
            position_[axis].push_back((i + 0.5) * spacing + offset);
        }
    }
}

/**
 * Bin phase space sample into spatial cells
 */
template <int dimension, typename float_type>
void binned_phase_space<dimension, float_type>::acquire()
{
    if (sample_->step == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return;
    }

    // trigger update of phase space sample
    on_acquire_();

    LOG_TRACE("acquire sample");
    scoped_timer_type timer(runtime_.acquire);

    binned_phase_space_wrapper<dimension> const& kernel_wrapper =
        binned_phase_space_wrapper<dimension>::kernel;

    // compute 1-dimensional cell indices for particle positions
    cuda::configure(dim_.grid, dim_.block);
    kernel_wrapper.compute_cell_index(
        sample_->position_data(), g_cell_index_
      , sample_->cell_length, sample_->cell_origin, sample_->nbin, npart_
    );

    // generate permutation
    cuda::configure(dim_.grid, dim_.block);
    kernel_wrapper.generate_sequence(g_cell_permutation_, npart_);
    sort_(g_cell_index_, g_cell_permutation_);

    // compute global cell offsets in sorted particle list
    assert(binned_phase_space_kernel::PLACEHOLDER == -1U);
    cuda::memset(g_cell_offset_, 0xFF);
    cuda::configure(dim_.grid, dim_.block);
    kernel_wrapper.find_cell_offset(g_cell_index_, g_cell_offset_, npart_);

    // assign particles to cells
    cuda::vector<int> g_ret(1);
    cuda::memset(g_ret, EXIT_SUCCESS);
    cuda::configure(sample_->dim.grid, sample_->dim.block);
    kernel_wrapper.assign_cells(
        g_ret, g_cell_index_, g_cell_offset_, g_cell_permutation_
      , sample_->g_cell, npart_, sample_->cell_size
    );
    vector<int> h_ret(1); // do not use page-locked memory here as allocation is very slow
    cuda::copy(g_ret, h_ret);
    if (h_ret.front() != EXIT_SUCCESS) {
        throw runtime_error("overcrowded placeholders in phase space binning");
    }

    sample_->step = clock_->step();
}

template <int dimension, typename float_type>
void binned_phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("binned_phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                class_<binned_phase_space, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                    ]
                    .def_readonly("runtime", &binned_phase_space::runtime_)
            ]

          , def("binned_phase_space", &make_shared<binned_phase_space
                , shared_ptr<sample_type>
                , shared_ptr<box_type const>
                , shared_ptr<clock_type const>
                , unsigned int
                , shared_ptr<logger_type>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_binned_phase_space(lua_State* L)
{
    binned_phase_space<3, float>::luaopen(L);
    binned_phase_space<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class binned_phase_space<3, float>;
template class binned_phase_space<2, float>;

} // namespace gpu
} // namespace observables
} // namespace halmd
