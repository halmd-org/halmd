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

#ifndef HALMD_MDSIM_HOST_FORCES_PAIR_TRUNC_HPP
#define HALMD_MDSIM_HOST_FORCES_PAIR_TRUNC_HPP

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <lua.hpp>
#include <string>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * template class for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_trunc
  : public mdsim::host::force<dimension, float_type>
{
public:
    typedef mdsim::host::force<dimension, float_type> _Base;
    typedef typename _Base::vector_type vector_type;
    typedef typename _Base::matrix_type matrix_type;
    typedef typename _Base::stress_tensor_type stress_tensor_type;

    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef host::forces::smooth<dimension, float_type> smooth_type;
    typedef utility::profiler profiler_type;

    boost::shared_ptr<potential_type> potential;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    boost::shared_ptr<smooth_type> smooth;

    inline static void luaopen(lua_State* L);

    inline pair_trunc(
        boost::shared_ptr<potential_type> potential
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      // FIXME , boost::shared_ptr<smooth_type> smooth
    );
    inline virtual void compute();
    inline void register_runtimes(profiler_type& profiler);

    //! return potential cutoffs
    virtual matrix_type const& cutoff()
    {
        return potential->r_cut();
    }

    //! return average potential energy per particle
    virtual double potential_energy()
    {
        return en_pot_;
    }

    //! potential part of stress tensor
    virtual stress_tensor_type stress_tensor_pot()
    {
        return stress_pot_;
    }

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG(
        compute_, std::string("computation of ") + potential_type::name() + " forces"
    );

protected:
    /** average potential energy per particle */
    double en_pot_;
    /** potential part of stress tensor */
    stress_tensor_type stress_pot_;

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
{}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * Compute pair forces, potential energy, and potential part of stress tensor
 */
template <int dimension, typename float_type, typename potential_type>
void pair_trunc<dimension, float_type, potential_type>::compute()
{
    scoped_timer<timer> timer_(boost::fusion::at_key<compute_>(runtime_));

    // initialise particle forces to zero
    std::fill(particle->f.begin(), particle->f.end(), 0);

    // initialise potential energy and stress tensor
    en_pot_ = 0;
    stress_pot_ = 0;

    for (size_t i = 0; i < particle->nbox; ++i) {
        // calculate pairwise Lennard-Jones force with neighbour particles
        BOOST_FOREACH(size_t j, particle->neighbour[i]) {
            // particle distance vector
            vector_type r = particle->r[i] - particle->r[j];
            box->reduce_periodic(r);
            // particle types
            unsigned a = particle->type[i];
            unsigned b = particle->type[j];
            // squared particle distance
            float_type rr = inner_prod(r, r);

            // truncate potential at cutoff length
            if (rr >= potential->rr_cut(a, b))
                continue;

            float_type fval, en_pot;
            boost::tie(fval, en_pot) = (*potential)(rr, a, b);

            // optionally smooth potential yielding continuous 2nd derivative
            // FIXME test performance of template versus runtime bool
            if (smooth) {
                smooth->compute(std::sqrt(rr), potential->r_cut(a, b), fval, en_pot);
            }

            // add force contribution to both particles
            particle->f[i] += r * fval;
            particle->f[j] -= r * fval;

            // add contribution to potential energy
            en_pot_ += en_pot;

            // ... and potential part of stress tensor
            stress_pot_ += fval * make_stress_tensor(rr, r);
        }
    }

    en_pot_ /= particle->nbox;
    stress_pot_ /= particle->nbox;

    // ensure that system is still in valid state
    if (isinf(en_pot_)) {
        throw std::runtime_error("Potential energy diverged");
    }
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
            namespace_("host")
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

}}} // namespace mdsim::host::forces

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_PAIR_TRUNC_HPP */
