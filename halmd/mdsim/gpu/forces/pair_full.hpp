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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/force.hpp>
#include <halmd/mdsim/gpu/forces/pair_full_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

#include <memory>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

/**
 * class template for modules implementing short ranged potential forces
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
    typedef typename _Base::stress_pot_type stress_pot_type;
    typedef typename _Base::stress_pot_array_type stress_pot_array_type;
    typedef typename _Base::hypervirial_array_type hypervirial_array_type;

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
    typedef typename potential_type::gpu_potential_type gpu_potential_type;
    typedef pair_full_wrapper<dimension, gpu_potential_type> gpu_wrapper;

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
    /** potential part of stress tensor for each particle  */
    cache<stress_pot_array_type> stress_pot_;
    /** hypervirial per particle */
    cache<hypervirial_array_type> hypervirial_;

    /** cache observer of net force per particle */
    cache<> net_force_cache_;
    /** cache observer of potential energy per particle */
    cache<> en_pot_cache_;
    /** cache observer of potential part of stress tensor per particle */
    cache<> stress_pot_cache_;
    /** cache observer of hypervirial per particle */
    cache<> hypervirial_cache_;

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
    auto net_force = make_cache_mutable(net_force_);
    auto en_pot = make_cache_mutable(en_pot_);
    auto stress_pot = make_cache_mutable(stress_pot_);
    auto hypervirial = make_cache_mutable(hypervirial_);

    net_force->reserve(particle_->dim.threads());
    en_pot->reserve(particle_->dim.threads());
    //
    // The GPU stores the stress tensor elements in column-major order to
    // optimise access patterns for coalescable access. Increase capacity of
    // GPU array such that there are 4 (6) in 2D (3D) elements per particle
    // available, although stress_pot_->size() still returns the number of
    // particles.
    //
    stress_pot->reserve(stress_pot_type::static_size * particle_->dim.threads());
    hypervirial->reserve(particle_->dim.threads());
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::net_force_array_type> const&
pair_full<dimension, float_type, potential_type>::net_force()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (net_force_cache_ != position_cache) {
        compute();
        net_force_cache_ = position_cache;
    }

    return net_force_;
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::en_pot_array_type> const&
pair_full<dimension, float_type, potential_type>::en_pot()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (en_pot_cache_ != position_cache) {
        compute_aux();
        net_force_cache_ = position_cache;
        en_pot_cache_ = position_cache;
        stress_pot_cache_ = position_cache;
        hypervirial_cache_ = position_cache;
    }

    return en_pot_;
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::stress_pot_array_type> const&
pair_full<dimension, float_type, potential_type>::stress_pot()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (stress_pot_cache_ != position_cache) {
        compute_aux();
        net_force_cache_ = position_cache;
        en_pot_cache_ = position_cache;
        stress_pot_cache_ = position_cache;
        hypervirial_cache_ = position_cache;
    }

    return stress_pot_;
}

template <int dimension, typename float_type, typename potential_type>
cache<typename pair_full<dimension, float_type, potential_type>::hypervirial_array_type> const&
pair_full<dimension, float_type, potential_type>::hypervirial()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (hypervirial_cache_ != position_cache) {
        compute_aux();
        net_force_cache_ = position_cache;
        en_pot_cache_ = position_cache;
        stress_pot_cache_ = position_cache;
        hypervirial_cache_ = position_cache;
    }

    return hypervirial_;
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::compute()
{
    LOG_TRACE("compute forces");

    position_array_type const& position = read_cache(particle_->position());
    auto net_force = make_cache_mutable(net_force_);

    scoped_timer_type timer(runtime_.compute);

    potential_->bind_textures();

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    gpu_wrapper::kernel.compute(
        &*net_force->begin()
        , &*position.begin()
        , nullptr
        , nullptr
        , nullptr
        , particle_->nparticle()
        , particle_->nspecies()
        , particle_->nspecies()
        , static_cast<position_type>(box_->length())
    );
    cuda::thread::synchronize();
}

template <int dimension, typename float_type, typename potential_type>
inline void pair_full<dimension, float_type, potential_type>::compute_aux()
{
    LOG_TRACE("compute forces with auxiliary variables");

    position_array_type const& position = read_cache(particle_->position());
    auto net_force = make_cache_mutable(net_force_);
    auto en_pot = make_cache_mutable(en_pot_);
    auto stress_pot = make_cache_mutable(stress_pot_);
    auto hypervirial = make_cache_mutable(hypervirial_);

    scoped_timer_type timer(runtime_.compute);

    potential_->bind_textures();

    cuda::configure(particle_->dim.grid, particle_->dim.block);
    gpu_wrapper::kernel.compute_aux(
        &*net_force->begin()
        , &*position.begin()
        , &*en_pot->begin()
        , &*stress_pot->begin()
        , &*hypervirial->begin()
        , particle_->nparticle()
        , particle_->nspecies()
        , particle_->nspecies()
        , static_cast<position_type>(box_->length())
    );
    cuda::thread::synchronize();
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
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_FORCES_PAIR_FULL_HPP */
