/*
 * Copyright © 2011-2012  Felix Höfling and Peter Colberg
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

#include <halmd/observables/gpu/density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension, typename float_type>
density_mode<dimension, float_type>::density_mode(
    shared_ptr<wavevector_type const> wavevector
  , shared_ptr<clock_type const> clock
  , shared_ptr<logger_type> logger
)
    // dependency injection
  : wavevector_(wavevector)
  , clock_(clock)
  , logger_(logger)
    // member initialisation
  , nq_(wavevector_->value().size())
  , dim_(50, 64 << DEVICE_SCALE) // at most 512 threads per block
    // memory allocation
  , g_q_(nq_)
  , g_sin_block_(nq_ * dim_.blocks_per_grid()), g_cos_block_(nq_ * dim_.blocks_per_grid())
  , g_sin_(nq_), g_cos_(nq_)
  , h_sin_(nq_), h_cos_(nq_)
{
    LOG_DEBUG(
        "[density_mode] CUDA configuration: " << dim_.blocks_per_grid() << " blocks of "
     << dim_.threads_per_block() << " threads each"
    );
    // copy wavevectors to CUDA device
    try {
        cuda::host::vector<gpu_vector_type> q(nq_);
        for (unsigned int i = 0; i < nq_; ++i) {
            // select wavevector from pair (wavenumber, wavevector),
            // cast from fixed_vector<double, ...> to fixed_vector<float, ...>
            // and finally to gpu_vector_type (float4 or float2)
            q[i] = static_cast<vector_type>(wavevector_->value()[i].second);
        }
        cuda::copy(q, g_q_);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to initialise device constants");
        throw;
    }

    // copy parameters to CUDA device
    try {
        cuda::copy(nq_, wrapper_type::kernel.nq);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to initialise device constants");
        throw;
    }
}

/**
 * Acquire sample of all density modes from phase space sample
 */
template <int dimension, typename float_type>
shared_ptr<typename density_mode<dimension, float_type>::sample_type const>
density_mode<dimension, float_type>::acquire(phase_space_type const& phase_space)
{
    scoped_timer_type timer(runtime_.acquire);

    if (rho_sample_ && rho_sample_->step == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return rho_sample_;
    }

    LOG_TRACE("acquire sample");

    if (phase_space.step != clock_->step()) {
        throw logic_error("GPU phase space sample was not updated");
    }

    // re-allocate memory which allows modules (e.g., dynamics::blocking_scheme)
    // to hold a previous copy of the sample
    rho_sample_ = make_shared<sample_type>(nq_);

    assert(rho_sample_->rho->size() == nq_);

    // compute density modes
    mode_vector_type& rho = *rho_sample_->rho; //< dereference shared_ptr
    try {
        cuda::configure(dim_.grid, dim_.block);
        wrapper_type::kernel.q.bind(g_q_);

        // compute exp(i q·r) for all wavevector/particle pairs and perform block sums
        wrapper_type::kernel.compute(
            *phase_space.r, phase_space.r->size()
          , g_sin_block_, g_cos_block_);
        cuda::thread::synchronize();

        // finalise block sums for each wavevector
        cuda::configure(
            nq_                        // #blocks: one per wavevector
          , dim_.block                 // #threads per block, must be a power of 2
        );
        wrapper_type::kernel.finalise(g_sin_block_, g_cos_block_, g_sin_, g_cos_, dim_.blocks_per_grid());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to compute density modes on GPU");
        throw;
    }

    // copy data from device and store in density_mode sample
    cuda::copy(g_sin_, h_sin_);
    cuda::copy(g_cos_, h_cos_);
    for (unsigned int i = 0; i < nq_; ++i) {
        rho[i] = mode_type(h_cos_[i], -h_sin_[i]);
    }
    rho_sample_->step = clock_->step();

    return rho_sample_;
}

template <typename sample_type, typename density_mode_type, typename slot_type>
static shared_ptr<sample_type const>
acquire(shared_ptr<density_mode_type> density_mode, slot_type const& phase_space)
{
    return density_mode->acquire(*phase_space());
}

template <typename sample_type, typename density_mode_type, typename slot_type>
static function<shared_ptr<sample_type const> ()>
wrap_acquire(shared_ptr<density_mode_type> density_mode, slot_type const& phase_space)
{
    return bind(&acquire<sample_type, density_mode_type, slot_type>, density_mode, phase_space);
}

template <int dimension, typename float_type>
void density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    typedef function<shared_ptr<phase_space_type const> ()> slot_type;

    using namespace luabind;
    static string class_name("density_mode_" + lexical_cast<string>(dimension));
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<density_mode>(class_name.c_str())
                .def("acquire", &wrap_acquire<sample_type, density_mode, slot_type>)
                .property("wavevector", &density_mode::wavevector)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("acquire", &runtime::acquire)
                ]
                .def_readonly("runtime", &density_mode::runtime_)

          , def("density_mode", &make_shared<density_mode
              , shared_ptr<wavevector_type const>
              , shared_ptr<clock_type const>
              , shared_ptr<logger_type>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_density_mode(lua_State* L)
{
    density_mode<3, float>::luaopen(L);
    density_mode<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class density_mode<3, float>;
template class density_mode<2, float>;

}}  // namespace observables::gpu

}  // namespace halmd
