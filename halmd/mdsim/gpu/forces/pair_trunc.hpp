/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2013      Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/forces/trunc/discontinuous.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.hpp>
#include <halmd/mdsim/gpu/neighbour.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

#include <memory>
#include <tuple>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing short ranged potential forces
 */
template <int dimension, typename float_type, typename potential_type, typename trunc_type = mdsim::forces::trunc::discontinuous>
class pair_trunc
  : public force<dimension, float_type>
{
private:
    typedef force<dimension, float_type> _Base;

public:
    typedef typename _Base::net_force_array_type net_force_array_type;
    typedef typename _Base::en_pot_array_type en_pot_array_type;
    typedef typename _Base::stress_pot_type stress_pot_type;
    typedef typename _Base::stress_pot_array_type stress_pot_array_type;
    typedef typename _Base::hypervirial_array_type hypervirial_array_type;

    typedef particle<dimension, float_type> particle_type;
    typedef box<dimension> box_type;
    typedef neighbour neighbour_type;
    typedef logger logger_type;

    pair_trunc(
        std::shared_ptr<potential_type const> potential
      , std::shared_ptr<particle_type const> particle1
      , std::shared_ptr<particle_type const> particle2
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<neighbour_type> neighbour
      , std::shared_ptr<trunc_type const> trunc = std::make_shared<trunc_type>()
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
    typedef typename neighbour_type::array_type neighbour_array_type;
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_trunc_wrapper<dimension, gpu_potential_type, trunc_type> gpu_wrapper;

    /** compute forces */
    void compute();
    /** compute forces with auxiliary variables */
    void compute_aux();

    /** pair potential */
    std::shared_ptr<potential_type const> potential_;
    /** state of first system */
    std::shared_ptr<particle_type const> particle1_;
    /** state of second system */
    std::shared_ptr<particle_type const> particle2_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** neighbour lists */
    std::shared_ptr<neighbour_type> neighbour_;
    /** smoothing functor */
    std::shared_ptr<trunc_type const> trunc_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** net force per particle */
    cache<net_force_array_type> net_force_;
    /** potential energy per particle */
    cache<en_pot_array_type> en_pot_;
    /** potential part of stress tensor for each particle */
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

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
pair_trunc<dimension, float_type, potential_type, trunc_type>::pair_trunc(
    std::shared_ptr<potential_type const> potential
  , std::shared_ptr<particle_type const> particle1
  , std::shared_ptr<particle_type const> particle2
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<neighbour_type> neighbour
  , std::shared_ptr<trunc_type const> trunc
  , std::shared_ptr<logger_type> logger
)
  : potential_(potential)
  , particle1_(particle1)
  , particle2_(particle2)
  , box_(box)
  , neighbour_(neighbour)
  , trunc_(trunc)
  , logger_(logger)
  , net_force_(particle1_->nparticle())
  , en_pot_(particle1_->nparticle())
  , stress_pot_(particle1_->nparticle())
  , hypervirial_(particle1_->nparticle())
{
    cache_proxy<net_force_array_type> net_force = net_force_;
    cache_proxy<en_pot_array_type> en_pot = en_pot_;
    cache_proxy<stress_pot_array_type> stress_pot = stress_pot_;
    cache_proxy<hypervirial_array_type> hypervirial = hypervirial_;

    net_force->reserve(particle1_->dim.threads());
    en_pot->reserve(particle1_->dim.threads());
    //
    // The GPU stores the stress tensor elements in column-major order to
    // optimise access patterns for coalescable access. Increase capacity of
    // GPU array such that there are 4 (6) in 2D (3D) elements per particle
    // available, although stress_pot_->size() still returns the number of
    // particles.
    //
    stress_pot->reserve(stress_pot_type::static_size * particle1_->dim.threads());
    hypervirial->reserve(particle1_->dim.threads());
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
cache<typename pair_trunc<dimension, float_type, potential_type, trunc_type>::net_force_array_type> const&
pair_trunc<dimension, float_type, potential_type, trunc_type>::net_force()
{
    cache<position_array_type> const& position1_cache = particle1_->position();
    cache<position_array_type> const& position2_cache = particle2_->position();

    if (net_force_cache_ != std::tie(position1_cache, position2_cache)) {
        compute();
        net_force_cache_ = std::tie(position1_cache, position2_cache);
    }
    return net_force_;
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
cache<typename pair_trunc<dimension, float_type, potential_type, trunc_type>::en_pot_array_type> const&
pair_trunc<dimension, float_type, potential_type, trunc_type>::en_pot()
{
    cache<position_array_type> const& position1_cache = particle1_->position();
    cache<position_array_type> const& position2_cache = particle2_->position();

    if (en_pot_cache_ != std::tie(position1_cache, position2_cache)) {
        compute_aux();
        net_force_cache_ = std::tie(position1_cache, position2_cache);
        en_pot_cache_ = net_force_cache_;
        stress_pot_cache_ = net_force_cache_;
        hypervirial_cache_ = net_force_cache_;
    }
    return en_pot_;
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
cache<typename pair_trunc<dimension, float_type, potential_type, trunc_type>::stress_pot_array_type> const&
pair_trunc<dimension, float_type, potential_type, trunc_type>::stress_pot()
{
    cache<position_array_type> const& position1_cache = particle1_->position();
    cache<position_array_type> const& position2_cache = particle2_->position();

    if (stress_pot_cache_ != std::tie(position1_cache, position2_cache)) {
        compute_aux();
        net_force_cache_ = std::tie(position1_cache, position2_cache);
        en_pot_cache_ = net_force_cache_;
        stress_pot_cache_ = net_force_cache_;
        hypervirial_cache_ = net_force_cache_;
    }
    return stress_pot_;
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
cache<typename pair_trunc<dimension, float_type, potential_type, trunc_type>::hypervirial_array_type> const&
pair_trunc<dimension, float_type, potential_type, trunc_type>::hypervirial()
{
    cache<position_array_type> const& position1_cache = particle1_->position();
    cache<position_array_type> const& position2_cache = particle2_->position();

    if (hypervirial_cache_ != std::tie(position1_cache, position2_cache)) {
        compute_aux();
        net_force_cache_ = std::tie(position1_cache, position2_cache);
        en_pot_cache_ = net_force_cache_;
        stress_pot_cache_ = net_force_cache_;
        hypervirial_cache_ = net_force_cache_;
    }
    return hypervirial_;
}


template <int dimension, typename float_type, typename potential_type, typename trunc_type>
inline void pair_trunc<dimension, float_type, potential_type, trunc_type>::compute()
{
    cache_proxy<position_array_type const> position1 = particle1_->position();
    cache_proxy<position_array_type const> position2 = particle2_->position();
    cache_proxy<neighbour_array_type const> g_neighbour = neighbour_->g_neighbour();
    cache_proxy<net_force_array_type> net_force = net_force_;

    LOG_TRACE("compute forces");

    scoped_timer_type timer(runtime_.compute);

    gpu_wrapper::kernel.r2.bind(*position2);
    potential_->bind_textures();

    cuda::configure(particle1_->dim.grid, particle1_->dim.block);
    gpu_wrapper::kernel.compute(
        &*position1->begin()
      , &*net_force->begin()
      , &*g_neighbour->begin()
      , neighbour_->size()
      , neighbour_->stride()
      , nullptr
      , nullptr
      , nullptr
      , particle1_->nspecies()
      , particle2_->nspecies()
      , static_cast<position_type>(box_->length())
      , *trunc_
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type, typename trunc_type>
inline void pair_trunc<dimension, float_type, potential_type, trunc_type>::compute_aux()
{
    cache_proxy<position_array_type const> position1 = particle1_->position();
    cache_proxy<position_array_type const> position2 = particle2_->position();
    cache_proxy<neighbour_array_type const> g_neighbour = neighbour_->g_neighbour();
    cache_proxy<net_force_array_type> net_force = net_force_;
    cache_proxy<en_pot_array_type> en_pot = en_pot_;
    cache_proxy<stress_pot_array_type> stress_pot = stress_pot_;
    cache_proxy<hypervirial_array_type> hypervirial = hypervirial_;

    LOG_TRACE("compute forces with auxiliary variables");

    scoped_timer_type timer(runtime_.compute);

    gpu_wrapper::kernel.r2.bind(*position2);
    potential_->bind_textures();

    cuda::configure(particle1_->dim.grid, particle1_->dim.block);
    gpu_wrapper::kernel.compute_aux(
        &*position1->begin()
      , &*net_force->begin()
      , &*g_neighbour->begin()
      , neighbour_->size()
      , neighbour_->stride()
      , &*en_pot->begin()
      , &*stress_pot->begin()
      , &*hypervirial->begin()
      , particle1_->nspecies()
      , particle2_->nspecies()
      , static_cast<position_type>(box_->length())
      , *trunc_
    );
    cuda::thread::synchronize();
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
                class_<pair_trunc, _Base>()
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("compute", &runtime::compute)
                    ]
                    .def_readonly("runtime", &pair_trunc::runtime_)

              , def("pair_trunc", &std::make_shared<pair_trunc,
                    std::shared_ptr<potential_type const>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<neighbour_type>
                  , std::shared_ptr<trunc_type const>
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_HPP */
