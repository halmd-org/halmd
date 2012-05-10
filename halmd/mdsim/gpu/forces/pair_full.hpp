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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <lua.hpp>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/forces/pair_full_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_full
  : public mdsim::gpu::force<dimension, float_type>
{
public:
    typedef mdsim::gpu::force<dimension, float_type> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_stress_tensor_type gpu_stress_tensor_type;
    typedef gpu::particle<dimension, float> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_full_wrapper<dimension, gpu_potential_type> gpu_wrapper;

    inline static void luaopen(lua_State* L);

    inline pair_full(
        boost::shared_ptr<potential_type> potential
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    inline virtual void compute();

    /**
     * enable computation of auxiliary variables
     *
     * The flag is reset by the next call to compute().
     */
    virtual void aux_enable()
    {
        LOG_TRACE("enable computation of auxiliary variables");
        aux_flag_ = true;
    }

    //! returns potential energies of particles
    virtual cuda::vector<float> const& potential_energy() const
    {
        assert_aux_valid();
        return g_en_pot_;
    }

    /** potential part of stress tensors of particles */
    virtual cuda::vector<gpu_stress_tensor_type> const& stress_tensor_pot() const
    {
        assert_aux_valid();
        return g_stress_pot_;
    }

    //! returns hyper virial of particles
    virtual cuda::vector<float> const& hypervirial() const
    {
        assert_aux_valid();
        return g_hypervirial_;
    }

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };

    void assert_aux_valid() const
    {
        if (!aux_valid_) {
            throw std::logic_error("Auxiliary variables were not enabled in force module.");
        }
    }

    boost::shared_ptr<potential_type> potential_;
    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<box_type> box_;

    /** flag for switching the computation of auxiliary variables in function compute() */
    bool aux_flag_;
    /** flag indicates that the auxiliary variables were updated by the last call to compute() */
    bool aux_valid_;
    /** potential energy for each particle */
    cuda::vector<float> g_en_pot_;
    /** potential part of stress tensor for each particle */
    cuda::vector<gpu_stress_tensor_type> g_stress_pot_;
    /** hyper virial for each particle */
    cuda::vector<float> g_hypervirial_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type>
pair_full<dimension, float_type, potential_type>::pair_full(
    boost::shared_ptr<potential_type> potential
  , boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type> box
)
  // dependency injection
  : potential_(potential)
  , particle_(particle)
  , box_(box)
  // member initalisation
  , aux_flag_(false)          //< disable auxiliary variables by default
  , aux_valid_(false)
  // memory allocation
  , g_en_pot_(particle_->dim.threads())
  , g_stress_pot_(particle_->dim.threads())
  , g_hypervirial_(particle_->dim.threads())
{
    cuda::copy(particle_->nbox, gpu_wrapper::kernel.npart);
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
    aux_valid_ = aux_flag_;
    if (!aux_flag_) {
        gpu_wrapper::kernel.compute(
            particle_->g_f, particle_->g_r, g_en_pot_, g_stress_pot_, g_hypervirial_
          , static_cast<vector_type>(box_->length())
        );
    }
    else {
        gpu_wrapper::kernel.compute_aux(
            particle_->g_f, particle_->g_r, g_en_pot_, g_stress_pot_, g_hypervirial_
          , static_cast<vector_type>(box_->length())
        );
        aux_flag_ = false;
    }
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
static char const* module_name_wrapper(pair_full<dimension, float_type, potential_type> const&)
{
    return potential_type::module_name();
}

template <int dimension, typename float_type, typename potential_type>
void pair_full<dimension, float_type, potential_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
    using namespace luabind;
    static std::string class_name("pair_full_" + boost::lexical_cast<std::string>(dimension) + "_");
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
                        class_<pair_full, boost::shared_ptr<_Base_Base>, bases<_Base_Base, _Base> >(potential_type::module_name())
                            .def(constructor<
                                boost::shared_ptr<potential_type>
                              , boost::shared_ptr<particle_type>
                              , boost::shared_ptr<box_type>
                            >())
                            .property("module_name", &module_name_wrapper<dimension, float_type, potential_type>)
                            .scope
                            [
                                class_<runtime>("runtime")
                                    .def_readonly("compute", &runtime::compute)
                            ]
                            .def_readonly("runtime", &pair_full::runtime_)
                    ]
                ]
            ]

          , namespace_("forces")
            [
                def("pair_full", &boost::make_shared<pair_full,
                    boost::shared_ptr<potential_type>
                  , boost::shared_ptr<particle_type>
                  , boost::shared_ptr<box_type>
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
