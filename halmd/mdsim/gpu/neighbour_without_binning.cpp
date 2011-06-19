/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <boost/bind.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/neighbour_without_binning.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * construct neighbour list module
 *
 * @param particle mdsim::gpu::particle instance
 * @param box mdsim::box instance
 * @param cutoff force cutoff radius
 * @param skin neighbour list skin
 * @param cell_occupancy desired average cell occupancy
 */
template <int dimension, typename float_type>
neighbour_without_binning<dimension, float_type>::neighbour_without_binning(
    shared_ptr<particle_type const> particle1
  , shared_ptr<particle_type const> particle2
  , shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , double cell_occupancy
)
  // dependency injection
  : particle1_(particle1)
  , particle2_(particle2)
  , box_(box)
  // allocate parameters
  , r_skin_(skin)
  , rr_cut_skin_(particle2_->ntype, particle2_->ntype)
  , g_rr_cut_skin_(rr_cut_skin_.data().size())
  , nu_cell_(cell_occupancy) // FIXME neighbour list occupancy
{
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < particle2_->ntype; ++i) {
        for (size_t j = i; j < particle2_->ntype; ++j) {
            rr_cut_skin_(i, j) = std::pow(r_cut(i, j) + r_skin_, 2);
            r_cut_max = max(r_cut(i, j), r_cut_max);
        }
    }
    cuda::copy(rr_cut_skin_.data(), g_rr_cut_skin_);

    // volume of n-dimensional sphere with neighbour list radius
    // volume of unit sphere: V_d = π^(d/2) / Γ(1+d/2), Γ(1) = 1, Γ(1/2) = √π
    float_type unit_sphere[5] = {0, 2, M_PI, 4 * M_PI / 3, M_PI * M_PI / 2 };
    assert(dimension <= 4);
    float_type neighbour_sphere = unit_sphere[dimension] * pow(r_cut_max + r_skin_, dimension);
    // number of placeholders per neighbour list
    size_ = static_cast<size_t>(ceil(neighbour_sphere * (box_->density() / nu_cell_)));
    // at least cell_size (or warp_size?) placeholders
    // FIXME what is a sensible lower bound?
    // size_ = max(size_, binning_->cell_size());
    // number of neighbour lists
    stride_ = particle1_->dim.threads();
    // allocate neighbour lists
    g_neighbour_.resize(stride_ * size_);

    LOG("neighbour list skin: " << r_skin_);
    LOG("number of placeholders per neighbour list: " << size_);
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void neighbour_without_binning<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.update, "update", "neighbour lists update");
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void neighbour_without_binning<dimension, float_type>::update()
{
    LOG_TRACE("update neighbour lists");

    scoped_timer<timer> timer_(runtime_.update);

    // mark neighbour list placeholders as virtual particles
    cuda::memset(g_neighbour_, 0xFF);
    // build neighbour lists
    cuda::vector<int> g_overflow(1);
    cuda::host::vector<int> h_overflow(1);
    cuda::memset(g_overflow, 0);
    get_neighbour_without_binning_kernel<dimension>().rr_cut_skin.bind(g_rr_cut_skin_);
    cuda::configure(
        particle1_->dim.grid
      , particle1_->dim.block
      , particle1_->dim.threads_per_block() * (sizeof(unsigned int) + sizeof(vector_type))
    );
    get_neighbour_without_binning_kernel<dimension>().update(
        particle1_->g_r
      , particle1_->nbox
      , particle2_->g_r
      , particle2_->nbox
      , static_cast<vector_type>(box_->length())
      , g_neighbour_
      , size_
      , stride_
      , g_overflow
    );
    cuda::copy(g_overflow, h_overflow);
    cuda::thread::synchronize();
    if (h_overflow.front() > 0) {
        LOG_ERROR("[neighbour] failed to bin " << h_overflow.front() << " particles");
        throw runtime_error("neighbour list occupancy too large");
    }
}

template <typename neighbour_type>
typename signal<void ()>::slot_function_type
wrap_update(shared_ptr<neighbour_type> neighbour)
{
    return bind(&neighbour_type::update, neighbour);
}

template <int dimension, typename float_type>
void neighbour_without_binning<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("neighbour_without_binning_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                class_<neighbour_without_binning, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def(constructor<
                        shared_ptr<particle_type const>
                      , shared_ptr<particle_type const>
                      , shared_ptr<box_type const>
                      , matrix_type const&
                      , double
                      , double
                    >())
                    .def("register_runtimes", &neighbour_without_binning::register_runtimes)
                    .property("update", &wrap_update<neighbour_without_binning>)
                    .property("r_skin", &neighbour_without_binning::r_skin)
                    .property("cell_occupancy", &neighbour_without_binning::cell_occupancy)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_neighbour_without_binning(lua_State* L)
{
    neighbour_without_binning<3, float>::luaopen(L);
    neighbour_without_binning<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class neighbour_without_binning<3, float>;
template class neighbour_without_binning<2, float>;

} // namespace gpu
} // namespace mdsim
} // namespace halmd
