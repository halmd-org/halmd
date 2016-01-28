/*
 * Copyright © 2016      Daniel Kirchner
 * Copyright © 2010-2012 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/config.hpp>

#include <halmd/algorithm/host/permute.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

#include <boost/function_output_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <exception>
#include <iterator>
#include <numeric>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
particle<dimension, float_type>::particle(size_type nparticle, unsigned int nspecies)
  : nparticle_(nparticle)
  , capacity_((nparticle + 128 - 1) & ~(128 - 1)) // round upwards to multiple of 128
  , nspecies_(std::max(nspecies, 1u))
  // enable auxiliary variables by default to allow sampling of initial state
  , force_zero_(true)
  , force_dirty_(true)
  , aux_dirty_(true)
  , aux_enabled_(true)
{
    // register and allocate named particle arrays
    auto position = make_cache_mutable(register_data<position_type>("position")->mutable_data());
    auto orientation = make_cache_mutable(register_data<orientation_type>("orientation")->mutable_data());
    auto image = make_cache_mutable(register_data<image_type>("image")->mutable_data());
    auto velocity = make_cache_mutable(register_data<velocity_type>("velocity")->mutable_data());
    auto id = make_cache_mutable(register_data<id_type>("id")->mutable_data());
    auto reverse_id = make_cache_mutable(register_data<reverse_id_type>("reverse_id")->mutable_data());
    auto species = make_cache_mutable(register_data<species_type>("species")->mutable_data());
    auto mass = make_cache_mutable(register_data<mass_type>("mass")->mutable_data());
    auto force = make_cache_mutable(register_data<force_type>("force", [this]() { this->update_force_(); })->mutable_data());
    auto en_pot = make_cache_mutable(register_data<en_pot_type>("potential_energy", [this]() { this->update_force_(true); })->mutable_data());
    auto stress_pot = make_cache_mutable(
            register_data<stress_pot_type>("potential_stress_tensor", [this]() { this->update_force_(true); })->mutable_data());

    // initialize particle arrays
    std::fill(position->begin(), position->end(), 0);
    std::fill(orientation->begin(), orientation->end(), 0);
    std::fill(image->begin(), image->end(), 0);
    std::fill(velocity->begin(), velocity->end(), 0);
    std::iota(id->begin(), id->begin() + nparticle_, 0);
    std::fill(id->begin() + nparticle_, id->end(), -1U);
    std::iota(reverse_id->begin(), reverse_id->begin() + nparticle_, 0);
    std::fill(reverse_id->begin() + nparticle_, reverse_id->end(), -1U);
    std::fill(species->begin(), species->begin() + nparticle_, 0);
    std::fill(species->begin() + nparticle_, species->end(), -1U);
    std::fill(mass->begin(), mass->end(), 1);
    std::fill(force->begin(), force->end(), 0);
    std::fill(en_pot->begin(), en_pot->end(), 0);
    std::fill(stress_pot->begin(), stress_pot->end(), 0);

    LOG("number of particles: " << nparticle_);
    LOG("number of particle species: " << nspecies_);
    LOG_DEBUG("capacity of data arrays: " << capacity_);
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::aux_enable()
{
    LOG_TRACE("enable computation of auxiliary variables");
    aux_enabled_ = true;
}

/**
 * Rearrange particles in memory according to an integer index sequence
 *
 * The neighbour lists must be rebuilt after calling this function!
 */
template <int dimension, typename float_type>
void particle<dimension, float_type>::rearrange(std::vector<unsigned int> const& index)
{
    scoped_timer_type timer(runtime_.rearrange);

    auto position = make_cache_mutable(mutable_data<position_type>("position"));
    auto orientation = make_cache_mutable(mutable_data<orientation_type>("orientation"));
    auto image = make_cache_mutable(mutable_data<image_type>("image"));
    auto velocity = make_cache_mutable(mutable_data<velocity_type>("velocity"));
    auto id = make_cache_mutable(mutable_data<id_type>("id"));
    auto reverse_id = make_cache_mutable(mutable_data<reverse_id_type>("reverse_id"));
    auto species = make_cache_mutable(mutable_data<species_type>("species"));
    auto mass = make_cache_mutable(mutable_data<mass_type>("mass"));

    permute(position->begin(), position->begin() + nparticle_, index.begin());
    permute(image->begin(), image->begin() + nparticle_, index.begin());
    permute(orientation->begin(), orientation->begin() + nparticle_, index.begin());
    permute(velocity->begin(), velocity->begin() + nparticle_, index.begin());
    permute(id->begin(), id->begin() + nparticle_, index.begin());
    permute(species->begin(), species->begin() + nparticle_, index.begin());
    permute(mass->begin(), mass->begin() + nparticle_, index.begin());
    // no permutation of forces

    // update reverse IDs
    for (unsigned int i = 0; i < nparticle_; ++i) {
        (*reverse_id)[(*id)[i]] = i;
    }
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::update_force_(bool with_aux)
{
    on_prepend_force_();          // ask force modules whether force/aux cache is dirty

    if (force_dirty_ || (with_aux && aux_dirty_)) {
        if (with_aux && aux_dirty_) {
            if (!force_dirty_) {
                LOG_WARNING_ONCE("auxiliary variables inactive in prior force computation, use aux_enable()");
            }
            aux_enabled_ = true;  // turn on computation of aux variables
        }
        LOG_TRACE("request force" << std::string(aux_enabled_ ? " and auxiliary variables" : ""));

        force_zero_ = true;       // tell first force module to reset the force
        on_force_();              // compute forces
        force_dirty_ = false;     // mark force cache as clean
        if (aux_enabled_) {
            aux_dirty_ = false;   // aux cache is clean only if requested
        }
        aux_enabled_ = false;     // disable aux variables for next call
    }
    on_append_force_();
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
}

template <typename particle_type>
static luaponte::object wrap_get(particle_type const& particle, lua_State* L, std::string const& name)
{
    return particle.get_array(name)->get_lua(L);
}

template <typename particle_type>
static void wrap_set(particle_type const& particle, std::string const& name, luaponte::object object)
{
    particle.get_array(name)->set_lua(object);
}

template <typename T>
static bool equal(std::shared_ptr<T const> self, std::shared_ptr<T const> other)
{
    // compare pointers of managed objects. I could not get working
    // owner_equal() or owner_before() with shared_ptr's passed from Lua
    return self == other;
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string class_name = "particle_" + std::to_string(dimension);
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<particle, std::shared_ptr<particle>>(class_name.c_str())
                    .def(constructor<size_type, unsigned int>())
                    .property("nparticle", &particle::nparticle)
                    .property("nspecies", &particle::nspecies)
                    .def("get", &wrap_get<particle>)
                    .def("set", &wrap_set<particle>)
                    .def("shift_velocity", &shift_velocity<particle>)
                    .def("shift_velocity_group", &shift_velocity_group<particle>)
                    .def("rescale_velocity", &rescale_velocity<particle>)
                    .def("rescale_velocity_group", &rescale_velocity_group<particle>)
                    .def("shift_rescale_velocity", &shift_rescale_velocity<particle>)
                    .def("shift_rescale_velocity_group", &shift_rescale_velocity_group<particle>)
                    .property("dimension", &wrap_dimension<dimension, float_type>)
                    .def("aux_enable", &particle::aux_enable)
                    .def("on_prepend_force", &particle::on_prepend_force)
                    .def("on_force", &particle::on_force)
                    .def("on_append_force", &particle::on_append_force)
                    .def("__eq", &equal<particle>) // operator= in Lua
                    .scope[
                        class_<runtime>("runtime")
                            .def_readonly("rearrange", &runtime::rearrange)
                    ]
                    .def_readonly("runtime", &particle::runtime_)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    particle<3, double>::luaopen(L);
    particle<2, double>::luaopen(L);
#else
    particle<3, float>::luaopen(L);
    particle<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class particle<3, double>;
template class particle<2, double>;
#else
template class particle<3, float>;
template class particle<2, float>;
#endif

} // namespace host
} // namespace mdsim
} // namespace halmd
