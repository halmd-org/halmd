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

#ifndef HALMD_MDSIM_HOST_FORCES_PAIR_FULL_HPP
#define HALMD_MDSIM_HOST_FORCES_PAIR_FULL_HPP

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <lua.hpp>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

/**
 * template class for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_full
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename potential_type::matrix_type matrix_type;

    typedef mdsim::box<dimension> box_type;
    typedef host::forces::smooth<dimension, float_type> smooth_type;

    static void luaopen(lua_State* L);

    pair_full(
        boost::shared_ptr<potential_type> potential
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      // FIXME , boost::shared_ptr<smooth_type> smooth
    );
    void compute();

    //! return potential cutoffs
    virtual matrix_type const& r_cut()
    {
        return potential_->r_cut();
    }

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };

    boost::shared_ptr<potential_type> potential_;
    boost::shared_ptr<particle_type> particle_;
    boost::shared_ptr<box_type> box_;
    boost::shared_ptr<smooth_type> smooth_;

    /** profiling runtime accumulators */
    runtime runtime_;

    template <bool compute_aux1, bool compute_aux2>
    inline void compute_aux();
};

template <int dimension, typename float_type, typename potential_type>
pair_full<dimension, float_type, potential_type>::pair_full(
    boost::shared_ptr<potential_type> potential
  , boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type> box
  // FIXME , boost::shared_ptr<smooth_type> smooth
)
  // dependency injection
  : potential_(potential)
  , particle_(particle)
  , box_(box)
{}

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

    if (particle_->aux_valid()) {
        compute_aux<true, true>();
    }
    else {
        compute_aux<false, false>();
    }
}

template <int dimension, typename float_type, typename potential_type>
template <bool compute_aux1, bool compute_aux2>
void pair_full<dimension, float_type, potential_type>::compute_aux()
{
    typename particle_type::force_array_type& force1             = particle_->force();
    typename particle_type::force_array_type& force2             = particle_->force();
    typename particle_type::en_pot_array_type& en_pot1           = particle_->en_pot();
    typename particle_type::en_pot_array_type& en_pot2           = particle_->en_pot();
    typename particle_type::stress_pot_array_type& stress_pot1   = particle_->stress_pot();
    typename particle_type::stress_pot_array_type& stress_pot2   = particle_->stress_pot();
    typename particle_type::hypervirial_array_type& hypervirial1 = particle_->hypervirial();
    typename particle_type::hypervirial_array_type& hypervirial2 = particle_->hypervirial();

    for (size_t i = 0; i < particle_->nparticle(); ++i) {
        // calculate untruncated pairwise Lennard-Jones force with all other particles
        for (size_t j = 0; j < particle_->nparticle(); ++j) {
            // skip self-interaction
            if (i == j ) {
                continue;
            }
            // particle distance vector
            vector_type r = particle_->r[i] - particle_->r[j];
            box_->reduce_periodic(r);
            // particle types
            unsigned a = particle_->type[i];
            unsigned b = particle_->type[j];
            // squared particle distance
            float_type rr = inner_prod(r, r);

            float_type fval, pot, hvir;
            boost::tie(fval, pot, hvir) = (*potential_)(rr, a, b);

            // optionally smooth potential yielding continuous 2nd derivative
            // FIXME test performance of template versus runtime bool
            if (smooth_) {
                smooth_->compute(std::sqrt(rr), potential_->r_cut(a, b), fval, pot);
            }

            // add force contribution to both particles
            force1[i] += r * fval;
            force2[j] -= r * fval;

            // contribution to potential energy
            typename particle_type::en_pot_type en_pot = 0.5 * pot;
            // potential part of stress tensor
            typename particle_type::stress_pot_type stress_pot = 0.5 * fval * make_stress_tensor(rr, r);
            // contribution to hypervirial
            typename particle_type::hypervirial_type hypervirial = 0.5 * hvir / (dimension * dimension);

            // add contributions for first particle
            if (compute_aux1) {
                en_pot1[i]      += en_pot;
                stress_pot1[i]  += stress_pot;
                hypervirial1[i] += hypervirial;
            }

            // add contributions for second particle
            if (compute_aux2) {
                en_pot2[j]      += en_pot;
                stress_pot2[j]  += stress_pot;
                hypervirial2[j] += hypervirial;
            }
        }
    }
}

template <int dimension, typename float_type, typename potential_type>
static char const* module_name_wrapper(pair_full<dimension, float_type, potential_type> const&)
{
    return potential_type::module_name();
}

template <typename force_type>
static typename signal<void ()>::slot_function_type
wrap_compute(boost::shared_ptr<force_type> force)
{
    return boost::bind(&force_type::compute, force);
}

template <int dimension, typename float_type, typename potential_type>
void pair_full<dimension, float_type, potential_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static std::string class_name("pair_full_" + boost::lexical_cast<std::string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("forces")
                [
                    namespace_(class_name.c_str())
                    [
                        class_<pair_full>(potential_type::module_name())
                            .property("r_cut", &pair_full::r_cut)
                            .property("module_name", &module_name_wrapper<dimension, float_type, potential_type>)
                            .property("compute", &wrap_compute<pair_full>)
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
} // namespace host
} // namespace forces
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_PAIR_FULL_HPP */
