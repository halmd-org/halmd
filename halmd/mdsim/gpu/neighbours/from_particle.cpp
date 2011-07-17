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

#include <halmd/mdsim/gpu/neighbours/from_particle.hpp>
#include <halmd/mdsim/gpu/neighbours/from_particle_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

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
from_particle<dimension, float_type>::from_particle(
    shared_ptr<particle_type const> particle1
  , shared_ptr<particle_type const> particle2
  , shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , shared_ptr<logger> logger
  , double cell_occupancy
)
  // dependency injection
  : particle1_(particle1)
  , particle2_(particle2)
  , box_(box)
  , logger_(logger)
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
void from_particle<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.update, "update", "neighbour lists update");
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void from_particle<dimension, float_type>::update()
{
    // Emit on_update signal, which may be connected e.g. to the binning
    // update slot. We don't call binning::update directly, since the order
    // of calls is setup at the Lua level, and it allows us to pass binning
    // as a const dependency.
    on_update_();

    LOG_TRACE("update neighbour lists");

    scoped_timer<timer> timer_(runtime_.update);

    // mark neighbour list placeholders as virtual particles
    cuda::memset(g_neighbour_, 0xFF);
    // build neighbour lists
    cuda::vector<int> g_overflow(1);
    cuda::host::vector<int> h_overflow(1);
    cuda::memset(g_overflow, 0);
    get_from_particle_kernel<dimension>().rr_cut_skin.bind(g_rr_cut_skin_);
    cuda::configure(
        particle1_->dim.grid
      , particle1_->dim.block
      , particle1_->dim.threads_per_block() * (sizeof(unsigned int) + sizeof(vector_type))
    );
    get_from_particle_kernel<dimension>().update(
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
        LOG_ERROR("failed to bin " << h_overflow.front() << " particles");
        throw runtime_error("neighbour list occupancy too large");
    }
}

template <int dimension, typename float_type>
float_type from_particle<dimension, float_type>::defaults::occupancy() {
    return 0.4;
}

template <int dimension, typename float_type>
void from_particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("from_particle_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("neighbours")
                [
                    class_<from_particle, shared_ptr<mdsim::neighbour>, mdsim::neighbour>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<particle_type const>
                          , shared_ptr<particle_type const>
                          , shared_ptr<box_type const>
                          , matrix_type const&
                          , double
                          , shared_ptr<logger_type>
                          , double
                        >())
                        .def("register_runtimes", &from_particle::register_runtimes)
                        .property("r_skin", &from_particle::r_skin)
                        .property("cell_occupancy", &from_particle::cell_occupancy)
                        .scope[
                            class_<defaults>("defaults")
                                .scope
                                [
                                    def("occupancy", &defaults::occupancy)
                                ]
                        ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_neighbours_from_particle(lua_State* L)
{
    from_particle<3, float>::luaopen(L);
    from_particle<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class from_particle<3, float>;
template class from_particle<2, float>;

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd
