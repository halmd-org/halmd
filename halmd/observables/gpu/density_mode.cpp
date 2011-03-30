/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/foreach.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/observables/gpu/density_mode.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace gpu
{

template <int dimension, typename float_type>
density_mode<dimension, float_type>::density_mode(
    shared_ptr<phase_space_type> phase_space
  , shared_ptr<wavevector_type> wavevector
)
    // dependency injection
  : phase_space_(phase_space)
  , wavevector_(wavevector)
    // member initialisation
  , nq_(wavevector_->value().size())
  , dim_(50, 64 << DEVICE_SCALE) // at most 512 threads per block
    // memory allocation
  , rho_sample_(phase_space_->sample->r.size(), nq_)
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
        LOG_ERROR("[density_mode] failed to initialise device constants");
        throw;
    }

    // copy parameters to CUDA device
    try {
        cuda::copy(nq_, wrapper_type::kernel.nq);
    }
    catch (cuda::error const&) {
        LOG_ERROR("[density_mode] failed to initialise device constants");
        throw;
    }
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void density_mode<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.sample, "sample", "computation of density modes");
}

/**
 * Acquire sample of all density modes from phase space sample
 */
template <int dimension, typename float_type>
void density_mode<dimension, float_type>::acquire(double time)
{
    scoped_timer<timer> timer_(runtime_.sample);

    if (rho_sample_.time == time) {
        LOG_TRACE("[density_mode] sample is up to date");
        return;
    }

    typedef typename phase_space_type::sample_type::sample_vector_ptr positions_vector_ptr_type;
    typedef typename density_mode_sample_type::mode_vector_type mode_vector_type;

    // trigger update of phase space sample
    on_acquire_(time);

    LOG_TRACE("[density_mode] acquire sample");

    if (phase_space_->sample->time != time) {
        throw logic_error("GPU phase space sample was not updated");
    }

    // compute density modes separately for each particle type
    // 1st loop: iterate over particle types
    unsigned int type = 0;
    BOOST_FOREACH (positions_vector_ptr_type const r_sample, phase_space_->sample->r) {
        mode_vector_type& rho = *rho_sample_.rho[type]; //< dereference shared_ptr
        if (type ==  0) {
            try {
                cuda::configure(dim_.grid, dim_.block);
                wrapper_type::kernel.q.bind(g_q_);

                // compute exp(i q·r) for all wavevector/particle pairs and perform block sums
                // FIXME pass r_sample->r[type] instead of particle->g_r
                wrapper_type::kernel.compute(
                    phase_space_->particle->g_r, phase_space_->particle->nbox
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
        }
        ++type;
    }
    rho_sample_.time = time;
}

template <int dimension, typename float_type>
void density_mode<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("density_mode_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("observables")
            [
                namespace_("gpu")
                [
                    class_<density_mode, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<phase_space_type>
                          , shared_ptr<wavevector_type>
                        >())
                        .def("register_runtimes", &density_mode::register_runtimes)
                ]
            ]
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
