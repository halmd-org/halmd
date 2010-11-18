/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
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
    typedef typename _Base::matrix_type matrix_type;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::gpu_stress_tensor_type gpu_stress_tensor_type;
    typedef gpu::particle<dimension, float> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef utility::profiler profiler_type;

    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_trunc_wrapper<dimension, gpu_potential_type> gpu_wrapper;

    boost::shared_ptr<potential_type> potential;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    inline static void luaopen(lua_State* L);

    inline pair_trunc(
        boost::shared_ptr<potential_type> potential
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    inline virtual void compute();
    inline void register_runtimes(profiler_type& profiler);

    //! returns potential cutoff distance
    virtual matrix_type const& cutoff()
    {
        return potential->r_cut();
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

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG(
        compute_, std::string("computation of ") + potential_type::name() + " forces"
    );

private:
    /** potential energy for each particle */
    cuda::vector<float> g_en_pot_;
    /** potential part of stress tensor for each particle */
    cuda::vector<gpu_stress_tensor_type> g_stress_pot_;

    boost::fusion::map<
        boost::fusion::pair<compute_, accumulator<double> >
    > runtime_;
};

template <int dimension, typename float_type, typename potential_type>
pair_trunc<dimension, float_type, potential_type>::pair_trunc(
    boost::shared_ptr<potential_type> potential
  , boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type> box
  // FIXME , boost::shared_ptr<smooth_type> smooth
)
  // dependency injection
  : potential(potential)
  , particle(particle)
  , box(box)
  // memory allocation
  , g_en_pot_(particle->dim.threads())
  , g_stress_pot_(particle->dim.threads())
{
    cuda::copy(static_cast<vector_type>(box->length()), gpu_wrapper::kernel.box_length);
}

/**
 * Compute pair forces, potential energy, and potential part of stress tensor
 */
template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::compute()
{
#ifdef USE_FORCE_DSFUN
#endif /* HALMD_VARIANT_FORCE_DSFUN */

    scoped_timer<timer> timer_(boost::fusion::at_key<compute_>(runtime_));

    cuda::copy(particle->neighbour_size, gpu_wrapper::kernel.neighbour_size);
    cuda::copy(particle->neighbour_stride, gpu_wrapper::kernel.neighbour_stride);
    gpu_wrapper::kernel.r.bind(particle->g_r);
    potential->bind_textures();

    cuda::configure(particle->dim.grid, particle->dim.block);
    gpu_wrapper::kernel.compute(particle->g_f, particle->g_neighbour, g_en_pot_, g_stress_pot_);
    cuda::thread::synchronize();
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
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
    module(L, "halmd_wrapper")
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
