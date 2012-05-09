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

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <lua.hpp>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/smoothers/nosmooth.hpp>
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
template <int dimension, typename float_type, typename potential_type, typename smooth_type = smoothers::nosmooth>
class pair_trunc
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef gpu::neighbour neighbour_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_trunc_wrapper<dimension, gpu_potential_type, smooth_type> gpu_wrapper;

    static void luaopen(lua_State* L);

    pair_trunc(
        boost::shared_ptr<potential_type> potential
      , boost::shared_ptr<particle_type> particle1
      , boost::shared_ptr<particle_type> particle2
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<neighbour_type const> neighbour
      , boost::shared_ptr<smooth_type const> smooth = boost::make_shared<smooth_type>()
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

    boost::shared_ptr<potential_type> potential_;
    boost::shared_ptr<particle_type> particle1_;
    boost::shared_ptr<particle_type> particle2_;
    boost::shared_ptr<box_type> box_;
    /** neighbour lists */
    boost::shared_ptr<neighbour_type const> neighbour_;
    /** smoothing functor */
    boost::shared_ptr<smooth_type const> smooth_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type, typename smooth_type>
pair_trunc<dimension, float_type, potential_type, smooth_type>::pair_trunc(
    boost::shared_ptr<potential_type> potential
  , boost::shared_ptr<particle_type> particle1
  , boost::shared_ptr<particle_type> particle2
  , boost::shared_ptr<box_type> box
  , boost::shared_ptr<neighbour_type const> neighbour
  , boost::shared_ptr<smooth_type const> smooth
)
  // dependency injection
  : potential_(potential)
  , particle1_(particle1)
  , particle2_(particle2)
  , box_(box)
  , neighbour_(neighbour)
  , smooth_(smooth) {}

/**
 * Compute pair forces and, if enabled, auxiliary variables,
 * i.e., potential energy, potential part of stress tensor
 *
 * Reset flag for auxiliary variables.
 */
template <int dimension, typename float_type, typename potential_type, typename smooth_type>
void pair_trunc<dimension, float_type, potential_type, smooth_type>::compute()
{
    scoped_timer_type timer(runtime_.compute);

    cuda::copy(neighbour_->size(), gpu_wrapper::kernel.neighbour_size);
    cuda::copy(neighbour_->stride(), gpu_wrapper::kernel.neighbour_stride);
    gpu_wrapper::kernel.r1.bind(particle1_->position());
    gpu_wrapper::kernel.r2.bind(particle2_->position());
    potential_->bind_textures();

    cuda::configure(particle1_->dim.grid, particle1_->dim.block);
    if (!particle1_->aux_valid()) {
        gpu_wrapper::kernel.compute(
            particle1_->force(), neighbour_->g_neighbour()
          , particle1_->en_pot(), particle1_->stress_pot(), particle1_->hypervirial()
          , particle1_->nspecies(), particle2_->nspecies()
          , static_cast<vector_type>(box_->length())
          , *smooth_
        );
    }
    else {
        gpu_wrapper::kernel.compute_aux(
            particle1_->force(), neighbour_->g_neighbour()
          , particle1_->en_pot(), particle1_->stress_pot(), particle1_->hypervirial()
          , particle1_->nspecies(), particle2_->nspecies()
          , static_cast<vector_type>(box_->length())
          , *smooth_
        );
    }
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type, typename smooth_type>
static char const* module_name_wrapper(pair_trunc<dimension, float_type, potential_type, smooth_type> const&)
{
    return potential_type::module_name();
}

template <typename force_type>
static typename signal<void ()>::slot_function_type
wrap_compute(boost::shared_ptr<force_type> force)
{
    return boost::bind(&force_type::compute, force);
}

template <int dimension, typename float_type, typename potential_type, typename smooth_type>
void pair_trunc<dimension, float_type, potential_type, smooth_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static std::string class_name("pair_trunc_" + boost::lexical_cast<std::string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("forces")
                [
                    namespace_(class_name.c_str())
                    [
                        class_<pair_trunc>(potential_type::module_name())
                            .property("module_name", &module_name_wrapper<dimension, float_type, potential_type, smooth_type>)
                            .property("compute", &wrap_compute<pair_trunc>)
                            .scope
                            [
                                class_<runtime>("runtime")
                                    .def_readonly("compute", &runtime::compute)
                            ]
                            .def_readonly("runtime", &pair_trunc::runtime_)
                    ]
                ]
            ]

          , namespace_("forces")
            [
                def("pair_trunc", &boost::make_shared<pair_trunc,
                    boost::shared_ptr<potential_type>
                  , boost::shared_ptr<particle_type>
                  , boost::shared_ptr<particle_type>
                  , boost::shared_ptr<box_type>
                  , boost::shared_ptr<neighbour_type const>
                  , boost::shared_ptr<smooth_type const>
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
