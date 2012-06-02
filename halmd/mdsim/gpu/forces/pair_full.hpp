/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP

#include <lua.hpp>
#include <memory>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/forces/pair_full_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_full
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_full_wrapper<dimension, gpu_potential_type> gpu_wrapper;

    static void luaopen(lua_State* L);

    pair_full(
        std::shared_ptr<potential_type> potential
      , std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type> box
    );
    void compute();

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };

    std::shared_ptr<potential_type> potential_;
    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<box_type> box_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type>
pair_full<dimension, float_type, potential_type>::pair_full(
    std::shared_ptr<potential_type> potential
  , std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type> box
)
  // dependency injection
  : potential_(potential)
  , particle_(particle)
  , box_(box)
{
    cuda::copy(particle_->nparticle(), gpu_wrapper::kernel.npart);
}

/**
 * Compute pair forces and, if enabled, auxiliary variables,
 * i.e., potential energy, potential part of stress tensor
 *
 * Reset flag for auxiliary variables.
 */
template <int dimension, typename float_type, typename potential_type>
void pair_full<dimension, float_type, potential_type>::compute()
{
    scoped_timer_type timer(runtime_.compute);

    potential_->bind_textures();

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    if (!particle_->aux_valid()) {
        gpu_wrapper::kernel.compute(
            particle_->force(), particle_->position()
          , particle_->en_pot(), particle_->stress_pot(), particle_->hypervirial()
          , particle_->nspecies(), particle_->nspecies()
          , static_cast<vector_type>(box_->length())
        );
    }
    else {
        gpu_wrapper::kernel.compute_aux(
            particle_->force(), particle_->position()
          , particle_->en_pot(), particle_->stress_pot(), particle_->hypervirial()
          , particle_->nspecies(), particle_->nspecies()
          , static_cast<vector_type>(box_->length())
        );
    }
    cuda::thread::synchronize();
}

template <typename force_type>
static std::function<void ()>
wrap_compute(std::shared_ptr<force_type> self)
{
    return [=]() {
        self->compute();
    };
}

template <int dimension, typename float_type, typename potential_type>
void pair_full<dimension, float_type, potential_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                class_<pair_full>()
                    .property("compute", &wrap_compute<pair_full>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &pair_full::runtime_)

              , def("pair_full", &std::make_shared<pair_full,
                    std::shared_ptr<potential_type>
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type>
                >)
            ]
        ]
    ];
}

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP */
