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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <lua.hpp>
#include <string>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

/**
 * class template for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_trunc
  : public mdsim::gpu::force<dimension, float_type>
{
public:
    typedef mdsim::gpu::force<dimension, float_type> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_stress_tensor_type gpu_stress_tensor_type;
    typedef gpu::particle<dimension, float> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef gpu::neighbour<dimension, float_type> neighbour_type;
    typedef utility::profiler profiler_type;

    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_trunc_wrapper<dimension, gpu_potential_type> gpu_wrapper;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type compute;
    };

    boost::shared_ptr<potential_type> potential;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    inline static void luaopen(lua_State* L);

    inline pair_trunc(
        boost::shared_ptr<potential_type> potential
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<neighbour_type const> neighbour
    );
    inline virtual void compute();
    inline void register_runtimes(profiler_type& profiler);

    //! enable computation of auxiliary variables
    virtual void aux_enable()
    {
        aux_flag_ = true;
    }

    //! disable computation of auxiliary variables
    virtual void aux_disable()
    {
        aux_flag_ = false;
    }

    //! return true if auxiliary variables are computed
    virtual bool aux_flag() const
    {
        return aux_flag_;
    }

    //! returns potential energies of particles
    virtual cuda::vector<float> const& potential_energy()
    {
        return g_en_pot_;
    }

    /** potential part of stress tensors of particles */
    virtual cuda::vector<gpu_stress_tensor_type> const& stress_tensor_pot()
    {
        return g_stress_pot_;
    }

    //! returns hyper virial of particles
    virtual cuda::vector<float> const& hypervirial()
    {
        return g_hypervirial_;
    }

private:
    /** neighbour lists */
    boost::shared_ptr<neighbour_type const> neighbour_;

    /** flag for switching the computation of auxiliary variables in function compute() */
    bool aux_flag_;
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
pair_trunc<dimension, float_type, potential_type>::pair_trunc(
    boost::shared_ptr<potential_type> potential
  , boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type> box
  , boost::shared_ptr<neighbour_type const> neighbour
  // FIXME , boost::shared_ptr<smooth_type> smooth
)
  // dependency injection
  : potential(potential)
  , particle(particle)
  , box(box)
  , neighbour_(neighbour)
  // member initalisation
  , aux_flag_(true)          //< enable everything by default
  // memory allocation
  , g_en_pot_(particle->dim.threads())
  , g_stress_pot_(particle->dim.threads())
  , g_hypervirial_(particle->dim.threads())
{
    cuda::copy(static_cast<vector_type>(box->length()), gpu_wrapper::kernel.box_length);
}

/**
 * Compute pair forces, potential energy, and potential part of stress tensor
 */
template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::compute()
{
    scoped_timer<timer> timer_(runtime_.compute);

    cuda::copy(neighbour_->size(), gpu_wrapper::kernel.neighbour_size);
    cuda::copy(neighbour_->stride(), gpu_wrapper::kernel.neighbour_stride);
    gpu_wrapper::kernel.r.bind(particle->g_r);
    potential->bind_textures();

    cuda::configure(particle->dim.grid, particle->dim.block);
    if (!aux_flag_) {
        gpu_wrapper::kernel.compute(
            particle->g_f, neighbour_->g_neighbour(), g_en_pot_, g_stress_pot_, g_hypervirial_
        );
    }
    else {
        gpu_wrapper::kernel.compute_aux(
            particle->g_f, neighbour_->g_neighbour(), g_en_pot_, g_stress_pot_, g_hypervirial_
        );
    }
    cuda::thread::synchronize();
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.compute, "compute", std::string("computation of ") + potential_type::name() + " forces");
}

template <int dimension, typename float_type, typename potential_type>
static char const* module_name_wrapper(pair_trunc<dimension, float_type, potential_type> const&)
{
    return potential_type::module_name();
}

template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base_Base;
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
                        class_<pair_trunc, boost::shared_ptr<_Base_Base>, bases<_Base_Base, _Base> >(potential_type::module_name())
                            .def(constructor<
                                boost::shared_ptr<potential_type>
                              , boost::shared_ptr<particle_type>
                              , boost::shared_ptr<box_type>
                              , boost::shared_ptr<neighbour_type const>
                            >())
                            .def("register_runtimes", &pair_trunc::register_runtimes)
                            .property("module_name", &module_name_wrapper<dimension, float_type, potential_type>)
                    ]
                ]
            ]
        ]
    ];
}

}}} // namespace mdsim::gpu::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_LJ_HPP */
