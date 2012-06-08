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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP

#include <lua.hpp>
#include <memory>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/forces/trunc/discontinuous.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
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
template <int dimension, typename float_type, typename potential_type, typename trunc_type = mdsim::forces::trunc::discontinuous>
class pair_trunc
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef gpu::neighbour neighbour_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_trunc_wrapper<dimension, gpu_potential_type, trunc_type> gpu_wrapper;

    static void luaopen(lua_State* L);

    pair_trunc(
        std::shared_ptr<potential_type> potential
      , std::shared_ptr<particle_type> particle1
      , std::shared_ptr<particle_type> particle2
      , std::shared_ptr<box_type> box
      , std::shared_ptr<neighbour_type const> neighbour
      , std::shared_ptr<trunc_type const> trunc = std::make_shared<trunc_type>()
    );
    void compute();

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::force_array_type force_array_type;
    typedef typename particle_type::en_pot_array_type en_pot_array_type;
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;
    typedef typename particle_type::hypervirial_array_type hypervirial_array_type;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };

    std::shared_ptr<potential_type> potential_;
    std::shared_ptr<particle_type> particle1_;
    std::shared_ptr<particle_type> particle2_;
    std::shared_ptr<box_type> box_;
    /** neighbour lists */
    std::shared_ptr<neighbour_type const> neighbour_;
    /** smoothing functor */
    std::shared_ptr<trunc_type const> trunc_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
pair_trunc<dimension, float_type, potential_type, trunc_type>::pair_trunc(
    std::shared_ptr<potential_type> potential
  , std::shared_ptr<particle_type> particle1
  , std::shared_ptr<particle_type> particle2
  , std::shared_ptr<box_type> box
  , std::shared_ptr<neighbour_type const> neighbour
  , std::shared_ptr<trunc_type const> trunc
)
  // dependency injection
  : potential_(potential)
  , particle1_(particle1)
  , particle2_(particle2)
  , box_(box)
  , neighbour_(neighbour)
  , trunc_(trunc) {}

/**
 * Compute pair forces and, if enabled, auxiliary variables,
 * i.e., potential energy, potential part of stress tensor
 *
 * Reset flag for auxiliary variables.
 */
template <int dimension, typename float_type, typename potential_type, typename trunc_type>
void pair_trunc<dimension, float_type, potential_type, trunc_type>::compute()
{
    cache_proxy<position_array_type const> position1 = particle1_->position();
    cache_proxy<position_array_type const> position2 = particle2_->position();
    cache_proxy<force_array_type> force1 = particle1_->force();
    cache_proxy<en_pot_array_type> en_pot1 = particle1_->en_pot();
    cache_proxy<stress_pot_array_type> stress_pot1 = particle1_->stress_pot();
    cache_proxy<hypervirial_array_type> hypervirial1 = particle1_->hypervirial();

    scoped_timer_type timer(runtime_.compute);

    cuda::copy(neighbour_->size(), gpu_wrapper::kernel.neighbour_size);
    cuda::copy(neighbour_->stride(), gpu_wrapper::kernel.neighbour_stride);
    gpu_wrapper::kernel.r1.bind(*position1);
    gpu_wrapper::kernel.r2.bind(*position2);
    potential_->bind_textures();

    cuda::configure(particle1_->dim.grid, particle1_->dim.block);
    if (!particle1_->aux_valid()) {
        gpu_wrapper::kernel.compute(
            &*force1->begin()
          , neighbour_->g_neighbour()
          , &*en_pot1->begin()
          , &*stress_pot1->begin()
          , &*hypervirial1->begin()
          , particle1_->nspecies()
          , particle2_->nspecies()
          , static_cast<vector_type>(box_->length())
          , *trunc_
        );
    }
    else {
        gpu_wrapper::kernel.compute_aux(
            &*force1->begin()
          , neighbour_->g_neighbour()
          , &*en_pot1->begin()
          , &*stress_pot1->begin()
          , &*hypervirial1->begin()
          , particle1_->nspecies()
          , particle2_->nspecies()
          , static_cast<vector_type>(box_->length())
          , *trunc_
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

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
void pair_trunc<dimension, float_type, potential_type, trunc_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                class_<pair_trunc>()
                    .property("compute", &wrap_compute<pair_trunc>)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &pair_trunc::runtime_)

              , def("pair_trunc", &std::make_shared<pair_trunc,
                    std::shared_ptr<potential_type>
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type>
                  , std::shared_ptr<neighbour_type const>
                  , std::shared_ptr<trunc_type const>
                >)
            ]
        ]
    ];
}

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LJ_HPP */
