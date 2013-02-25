/*
 * Copyright © 2010-2011 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

#include <memory>
#include <tuple>

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

/**
 * template class for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type>
class pair_full
  : public force<dimension, float_type>
{
private:
    typedef force<dimension, float_type> _Base;

public:
    typedef typename _Base::net_force_array_type net_force_array_type;
    typedef typename _Base::en_pot_array_type en_pot_array_type;
    typedef typename _Base::stress_pot_array_type stress_pot_array_type;
    typedef typename _Base::hypervirial_array_type hypervirial_array_type;
    typedef typename _Base::stress_pot_type stress_pot_type;

    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;
    typedef logger logger_type;

    pair_full(
        std::shared_ptr<potential_type const> potential
      , std::shared_ptr<particle_type const> particle
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /**
     * Returns const reference to net force per particle.
     */
    virtual cache<net_force_array_type> const& net_force();

    /**
     * Returns const reference to potential energy per particle.
     */
    virtual cache<en_pot_array_type> const& en_pot();

    /**
     * Returns const reference to potential part of stress tensor per particle.
     */
    virtual cache<stress_pot_array_type> const& stress_pot();

    /**
     * Returns const reference to hypervirial per particle.
     */
    virtual cache<hypervirial_array_type> const& hypervirial();

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::position_type position_type;
    typedef typename particle_type::species_array_type species_array_type;
    typedef typename particle_type::species_type species_type;
    typedef typename particle_type::size_type size_type;
    typedef typename _Base::en_pot_type en_pot_type;
    typedef typename _Base::hypervirial_type hypervirial_type;

    /** compute forces */
    void compute();
    /** compute forces with auxiliary variables */
    void compute_aux();

    /** pair potential */
    std::shared_ptr<potential_type const> potential_;
    /** system state */
    std::shared_ptr<particle_type const> particle_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** net force per particle */
    cache<net_force_array_type> net_force_;
    /** potential energy per particle */
    cache<en_pot_array_type> en_pot_;
    /** potential part of stress tensor of each particle */
    cache<stress_pot_array_type> stress_pot_;
    /** hypervirial per particle */
    cache<hypervirial_array_type> hypervirial_;

    /** cache observer of net force per particle */
    std::tuple<cache<>, cache<>> net_force_cache_;
    /** cache observer of potential energy per particle */
    std::tuple<cache<>, cache<>> en_pot_cache_;
    /** cache observer of potential part of stress tensor per particle */
    std::tuple<cache<>, cache<>> stress_pot_cache_;
    /** cache observer of hypervirial per particle */
    std::tuple<cache<>, cache<>> hypervirial_cache_;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type compute;
    };
    /** profiling runtime accumulators */
    runtime runtime_;
};

template <int dimension, typename float_type, typename potential_type>
pair_full<dimension, float_type, potential_type>::pair_full(
    std::shared_ptr<potential_type const> potential
  , std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger_type> logger
)
  : potential_(potential)
  , particle_(particle)
  , box_(box)
  , logger_(logger)
  , net_force_(particle_->nparticle())
  , en_pot_(particle_->nparticle())
  , stress_pot_(particle_->nparticle())
  , hypervirial_(particle_->nparticle())
{
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::net_force_array_type> const&
pair_full<dimension, float_type, potential_type>::net_force()
{
    cache<position_array_type> const& position_cache = particle_->position();
    cache<species_array_type> const& species_cache = particle_->species();

    if (net_force_cache_ != std::tie(position_cache, species_cache)) {
        LOG_TRACE("compute net force per particle");

        compute();

        net_force_cache_ = std::tie(position_cache, species_cache);
    }

    return net_force_;
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::en_pot_array_type> const&
pair_full<dimension, float_type, potential_type>::en_pot()
{
    cache<position_array_type> const& position_cache = particle_->position();
    cache<species_array_type> const& species_cache = particle_->species();

    if (en_pot_cache_ != std::tie(position_cache, species_cache)) {
        LOG_TRACE("compute potential energy per particle");

        compute_aux();

        net_force_cache_ = std::tie(position_cache, species_cache);
        en_pot_cache_ = net_force_cache_;
        stress_pot_cache_ = net_force_cache_;
        hypervirial_cache_ = net_force_cache_;
    }

    return en_pot_;
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::stress_pot_array_type> const&
pair_full<dimension, float_type, potential_type>::stress_pot()
{
    cache<position_array_type> const& position_cache = particle_->position();
    cache<species_array_type> const& species_cache = particle_->species();

    if (stress_pot_cache_ != std::tie(position_cache, species_cache)) {
        LOG_TRACE("compute potential part of stress tensor per particle");

        compute_aux();

        net_force_cache_ = std::tie(position_cache, species_cache);
        en_pot_cache_ = net_force_cache_;
        stress_pot_cache_ = net_force_cache_;
        hypervirial_cache_ = net_force_cache_;
    }

    return stress_pot_;
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::hypervirial_array_type> const&
pair_full<dimension, float_type, potential_type>::hypervirial()
{
    cache<position_array_type> const& position_cache = particle_->position();
    cache<species_array_type> const& species_cache = particle_->species();

    if (hypervirial_cache_ != std::tie(position_cache, species_cache)) {
        LOG_TRACE("compute hypervirial per particle");

        compute_aux();

        net_force_cache_ = std::tie(position_cache, species_cache);
        en_pot_cache_ = net_force_cache_;
        stress_pot_cache_ = net_force_cache_;
        hypervirial_cache_ = net_force_cache_;
    }

    return hypervirial_;
}


template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::compute()
{
    auto net_force = make_cache_mutable(net_force_);

    position_array_type const& position1 = read_cache(particle_->position());
    position_array_type const& position2 = read_cache(particle_->position());
    species_array_type const& species1   = *particle_->species();
    species_array_type const& species2   = *particle_->species();
    size_type nparticle1 = particle_->nparticle();
    size_type nparticle2 = particle_->nparticle();

    LOG_TRACE("compute forces");

    scoped_timer_type timer(runtime_.compute);

    std::fill(net_force->begin(), net_force->end(), 0);

    for (size_type i = 0; i < nparticle1; ++i) {
        // calculate untruncated pairwise Lennard-Jones force with all other particles
        for (size_type j = 0; j < nparticle2; ++j) {
            // skip self-interaction
            if (i == j ) {
                continue;
            }
            // particle distance vector
            position_type r = position1[i] - position2[j];
            box_->reduce_periodic(r);
            // particle types
            species_type a = species1[i];
            species_type b = species2[j];
            // squared particle distance
            float_type rr = inner_prod(r, r);

            float_type fval, pot, hvir;
            boost::tie(fval, pot, hvir) = (*potential_)(rr, a, b);

            // add force contribution to both particles
            (*net_force)[i] += r * fval;
            (*net_force)[j] -= r * fval;
        }
    }
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::compute_aux()
{
    cache_proxy<net_force_array_type> net_force     = net_force_;
    cache_proxy<en_pot_array_type> en_pot           = en_pot_;
    cache_proxy<stress_pot_array_type> stress_pot   = stress_pot_;
    auto hypervirial = make_cache_mutable(hypervirial_);

    position_array_type const& position1 = read_cache(particle_->position());
    position_array_type const& position2 = read_cache(particle_->position());
    species_array_type const& species1   = *particle_->species();
    species_array_type const& species2   = *particle_->species();
    size_type nparticle1 = particle_->nparticle();
    size_type nparticle2 = particle_->nparticle();

    LOG_TRACE("compute forces with auxiliary variables");

    scoped_timer_type timer(runtime_.compute);

    std::fill(net_force->begin(), net_force->end(), 0);
    std::fill(en_pot->begin(), en_pot->end(), 0);
    std::fill(stress_pot->begin(), stress_pot->end(), 0);
    std::fill(hypervirial->begin(), hypervirial->end(), 0);

    for (size_type i = 0; i < nparticle1; ++i) {
        // calculate untruncated pairwise Lennard-Jones force with all other particles
        for (size_type j = 0; j < nparticle2; ++j) {
            // skip self-interaction
            if (i == j ) {
                continue;
            }
            // particle distance vector
            position_type r = position1[i] - position2[j];
            box_->reduce_periodic(r);
            // particle types
            species_type a = species1[i];
            species_type b = species2[j];
            // squared particle distance
            float_type rr = inner_prod(r, r);

            float_type fval, pot, hvir;
            boost::tie(fval, pot, hvir) = (*potential_)(rr, a, b);

            // add force contribution to both particles
            (*net_force)[i] += r * fval;
            (*net_force)[j] -= r * fval;

            // contribution to potential energy
            en_pot_type en = 0.5 * pot;
            // potential part of stress tensor
            stress_pot_type stress = 0.5 * fval * make_stress_tensor(r);
            // contribution to hypervirial
            hypervirial_type hyper = 0.5 * hvir / (dimension * dimension);

            // add contributions for first particle
            (*en_pot)[i]      += en;
            (*stress_pot)[i]  += stress;
            (*hypervirial)[i] += hyper;

            // add contributions for second particle
            (*en_pot)[j]      += en;
            (*stress_pot)[j]  += stress;
            (*hypervirial)[j] += hyper;
        }
    }
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
                class_<pair_full, _Base>()
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &pair_full::runtime_)

              , def("pair_full", &std::make_shared<pair_full,
                    std::shared_ptr<potential_type const>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_FORCES_PAIR_FULL_HPP */
