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

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/power_law.hpp>
#include <halmd/numeric/pow.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

template <int dimension, typename float_type>
int const power_law<dimension, float_type>::default_index = 12;

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void power_law<dimension, float_type>::options(po::options_description& desc)
{
    desc.add_options()
        ("power-law-index", po::value<int>()->default_value(default_index),
         "index of soft power-law potential")
        ;
}

/**
 * Register option value types with Lua
 */
static __attribute__((constructor)) void register_option_converters()
{
    register_any_converter<int>();
}

/**
 * Initialize Lennard-Jones potential parameters
 */
template <int dimension, typename float_type>
power_law<dimension, float_type>::power_law(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
// FIXME  , shared_ptr<smooth_type> smooth
  , int index
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
)
  // dependency injection
  : particle(particle)
  , box(box)
// FIXME  , smooth(smooth)
  // allocate potential parameters
  , index_(index)
  , epsilon_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
  , sigma_(scalar_matrix<float_type>(particle->ntype, particle->ntype, 1))
  , r_cut_sigma_(particle->ntype, particle->ntype)
  , r_cut_(particle->ntype, particle->ntype)
  , rr_cut_(particle->ntype, particle->ntype)
  , sigma2_(particle->ntype, particle->ntype)
  , en_cut_(particle->ntype, particle->ntype)
{
    // FIXME support any number of types
    for (size_t i = 0; i < std::min(particle->ntype, 2U); ++i) {
        for (size_t j = i; j < std::min(particle->ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            r_cut_sigma_(i, j) = cutoff[i + j];
        }
    }

    // precalculate derived parameters
    for (size_t i = 0; i < particle->ntype; ++i) {
        for (size_t j = i; j < particle->ntype; ++j) {
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            // energy shift due to truncation at cutoff length
            en_cut_(i, j) = epsilon_(i, j) * std::pow(1 / r_cut_sigma_(i, j), static_cast<int>(index_));
        }
    }

    LOG("potential power-law index: n = " << index_);
    LOG("potential strength: ε = " << epsilon_);
    LOG("potential pair separation: σ = " << sigma_);
    LOG("potential cutoff length: r = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);
}

/**
 * Compute power-law forces
 *
 * call index-dependent template implementations
 * for efficiency of fixed_pow() function
 */
template <int dimension, typename float_type>
void power_law<dimension, float_type>::compute()
{
    switch (index_) {
        case 6:  compute_impl<6>();  break;
        case 12: compute_impl<12>(); break;
        case 24: compute_impl<24>(); break;
        case 48: compute_impl<48>(); break;
        default:
            LOG_WARNING_ONCE("Using non-optimised force routine for index " << index_);
            compute_impl<0>();
            break;
    }
}

template <int dimension, typename float_type>
template <int index>
void power_law<dimension, float_type>::compute_impl()
{
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
            size_t a = particle->type[i];
            size_t b = particle->type[j];
            // squared particle distance
            float_type rr = inner_prod(r, r);

            // truncate potential at cutoff length
            if (rr >= rr_cut_(a, b))
                continue;

            // compute power-law force in reduced units,
            // choose arbitrary index_ if template parameter index = 0
            float_type rni;
            if (index > 0) {
                rni = fixed_pow<index>(sigma_(a, b) / std::sqrt(rr));
            }
            else {
                rni = std::pow(sigma_(a, b) / std::sqrt(rr), index_);
            }
            float_type en_pot = epsilon_(a, b) * rni;      // U(r)
            float_type fval = (index > 0 ? index : index_) * en_pot / rr;
                                                           // F(r) / r
            en_pot -= en_cut_(a, b);                       // shift potential

            // optionally smooth potential yielding continuous 2nd derivative
            // FIXME test performance of template versus runtime bool
            // if (smooth) {
            //    smooth->compute(std::sqrt(rr), r_cut_(a, b), fval, en_pot);
            // }

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
    if (std::isinf(en_pot_)) {
        throw runtime_error("Potential energy diverged");
    }
}

template <typename T>
static void register_lua(lua_State* L, char const* class_name)
{
    typedef typename T::_Base _Base;
    typedef typename _Base::_Base _Base_Base;
    typedef typename T::particle_type particle_type;
    typedef typename T::box_type box_type;

    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("forces")
                    [
                        class_<T, shared_ptr<_Base_Base>, bases<_Base, _Base_Base> >(class_name)
                            .def(constructor<
                                shared_ptr<particle_type>
                              , shared_ptr<box_type>
                              , int
                              , array<float, 3> const&
                              , array<float, 3> const&
                              , array<float, 3> const&
                            >())
                            .scope
                            [
                                def("options", &T::options)
                            ]
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        bind(&register_lua<power_law<3, double> >, _1, "power_law_3_")
    ]
    [
        bind(&register_lua<power_law<2, double> >, _1, "power_law_2_")
    ];
#else
    [
        bind(&register_lua<power_law<3, float> >, _1, "power_law_3_")
    ]
    [
        bind(&register_lua<power_law<2, float> >, _1, "power_law_2_")
    ];
#endif
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class power_law<3, double>;
template class power_law<2, double>;
#else
template class power_law<3, float>;
template class power_law<2, float>;
#endif

}}} // namespace mdsim::host::forces

} // namespace halmd
