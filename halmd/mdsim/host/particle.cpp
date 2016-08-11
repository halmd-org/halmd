/*
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
  // allocate particle storage
  : nparticle_(nparticle)
  , nspecies_(std::max(nspecies, 1u))
  , position_(nparticle)
  , image_(nparticle)
  , velocity_(nparticle)
  , tag_(nparticle)
  , reverse_tag_(nparticle)
  , species_(nparticle)
  , mass_(nparticle)
  , force_(nparticle)
  , en_pot_(nparticle)
  , stress_pot_(nparticle)
  // enable auxiliary variables by default to allow sampling of initial state
  , force_zero_(true)
  , force_dirty_(true)
  , aux_dirty_(true)
  , aux_enabled_(true)
{
    auto position = make_cache_mutable(position_);
    auto image = make_cache_mutable(image_);
    auto velocity = make_cache_mutable(velocity_);
    auto tag = make_cache_mutable(tag_);
    auto reverse_tag = make_cache_mutable(reverse_tag_);
    auto species = make_cache_mutable(species_);
    auto mass = make_cache_mutable(mass_);
    auto force = make_cache_mutable(force_);
    auto en_pot = make_cache_mutable(en_pot_);
    auto stress_pot = make_cache_mutable(stress_pot_);

    std::fill(position->begin(), position->end(), 0);
    std::fill(image->begin(), image->end(), 0);
    std::fill(velocity->begin(), velocity->end(), 0);
    std::iota(tag->begin(), tag->end(), 0);
    std::iota(reverse_tag->begin(), reverse_tag->end(), 0);
    std::fill(species->begin(), species->end(), 0);
    std::fill(mass->begin(), mass->end(), 1);
    std::fill(force->begin(), force->end(), 0);
    std::fill(en_pot->begin(), en_pot->end(), 0);
    std::fill(stress_pot->begin(), stress_pot->end(), 0);

    LOG("number of particles: " << nparticle_);
    LOG("number of particle species: " << nspecies_);
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

    auto position = make_cache_mutable(position_);
    auto image = make_cache_mutable(image_);
    auto velocity = make_cache_mutable(velocity_);
    auto tag = make_cache_mutable(tag_);
    auto reverse_tag = make_cache_mutable(reverse_tag_);
    auto species = make_cache_mutable(species_);
    auto mass = make_cache_mutable(mass_);

    permute(position->begin(), position->end(), index.begin());
    permute(image->begin(), image->end(), index.begin());
    permute(velocity->begin(), velocity->end(), index.begin());
    permute(tag->begin(), tag->end(), index.begin());
    permute(species->begin(), species->end(), index.begin());
    permute(mass->begin(), mass->end(), index.begin());
    // no permutation of forces

    // update reverse tags
    for (unsigned int i = 0; i < nparticle_; ++i) {
        (*reverse_tag)[(*tag)[i]] = i;
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

template <typename particle_type>
static std::vector<typename particle_type::position_type>
wrap_get_position(particle_type const& self)
{
    std::vector<typename particle_type::position_type> output;
    output.reserve(self.nparticle());
    get_position(self, back_inserter(output));
    return std::move(output);
}

template <typename particle_type>
static void
wrap_set_position(particle_type& self, std::vector<typename particle_type::position_type> const& input)
{
    if (input.size() != self.nparticle()) {
        throw std::invalid_argument("input array size not equal to number of particles");
    }
    set_position(self, input.begin());
}

template <typename particle_type>
static std::vector<typename particle_type::image_type>
wrap_get_image(particle_type const& self)
{
    std::vector<typename particle_type::image_type> output;
    output.reserve(self.nparticle());
    get_image(self, back_inserter(output));
    return std::move(output);
}

template <typename particle_type>
static void
wrap_set_image(particle_type& self, std::vector<typename particle_type::image_type> const& input)
{
    if (input.size() != self.nparticle()) {
        throw std::invalid_argument("input array size not equal to number of particles");
    }
    set_image(self, input.begin());
}

template <typename particle_type>
static std::vector<typename particle_type::velocity_type>
wrap_get_velocity(particle_type const& self)
{
    std::vector<typename particle_type::velocity_type> output;
    output.reserve(self.nparticle());
    get_velocity(self, back_inserter(output));
    return std::move(output);
}

template <typename particle_type>
static void
wrap_set_velocity(particle_type& self, std::vector<typename particle_type::velocity_type> const& input)
{
    if (input.size() != self.nparticle()) {
        throw std::invalid_argument("input array size not equal to number of particles");
    }
    set_velocity(self, input.begin());
}

template <typename particle_type>
static std::vector<typename particle_type::tag_type>
wrap_get_tag(particle_type const& self)
{
    std::vector<typename particle_type::tag_type> output;
    output.reserve(self.nparticle());
    get_tag(self, back_inserter(output));
    return std::move(output);
}

template <typename particle_type>
static void
wrap_set_tag(particle_type& self, std::vector<typename particle_type::tag_type> const& input)
{
    typedef typename particle_type::tag_type tag_type;
    if (input.size() != self.nparticle()) {
        throw std::invalid_argument("input array size not equal to number of particles");
    }
    tag_type nparticle = self.nparticle();
    set_tag(
        self
      , boost::make_transform_iterator(input.begin(), [&](tag_type t) -> tag_type {
            if (t >= nparticle) {
                throw std::invalid_argument("invalid particle tag");
            }
            return t;
        })
    );
}

template <typename particle_type>
static std::vector<typename particle_type::reverse_tag_type>
wrap_get_reverse_tag(particle_type const& self)
{
    std::vector<typename particle_type::reverse_tag_type> output;
    output.reserve(self.nparticle());
    get_reverse_tag(self, back_inserter(output));
    return std::move(output);
}

template <typename particle_type>
static void
wrap_set_reverse_tag(particle_type& self, std::vector<typename particle_type::reverse_tag_type> const& input)
{
    typedef typename particle_type::reverse_tag_type reverse_tag_type;
    if (input.size() != self.nparticle()) {
        throw std::invalid_argument("input array size not equal to number of particles");
    }
    reverse_tag_type nparticle = self.nparticle();
    set_reverse_tag(
        self
      , boost::make_transform_iterator(input.begin(), [&](reverse_tag_type i) -> reverse_tag_type {
            if (i >= nparticle) {
                throw std::invalid_argument("invalid particle reverse tag");
            }
            return i;
        })
    );
}

template <typename particle_type>
static std::vector<typename particle_type::species_type>
wrap_get_species(particle_type const& self)
{
    std::vector<typename particle_type::species_type> output;
    output.reserve(self.nparticle());
    get_species(self, back_inserter(output));
    return std::move(output);
}

template <typename particle_type>
static void
wrap_set_species(particle_type& self, std::vector<typename particle_type::species_type> const& input)
{
    typedef typename particle_type::species_type species_type;
    if (input.size() != self.nparticle()) {
        throw std::invalid_argument("input array size not equal to number of particles");
    }
    species_type nspecies = self.nspecies();
    set_species(
        self
      , boost::make_transform_iterator(input.begin(), [&](species_type s) -> species_type {
            if (s >= nspecies) {
                throw std::invalid_argument("invalid particle species");
            }
            return s;
        })
    );
}

template <typename particle_type>
static std::vector<typename particle_type::mass_type>
wrap_get_mass(particle_type const& self)
{
    std::vector<typename particle_type::mass_type> output;
    output.reserve(self.nparticle());
    get_mass(self, back_inserter(output));
    return std::move(output);
}

template <typename particle_type>
static void
wrap_set_mass(particle_type& self, std::vector<typename particle_type::mass_type> const& input)
{
    if (input.size() != self.nparticle()) {
        throw std::invalid_argument("input array size not equal to number of particles");
    }
    set_mass(self, input.begin());
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::force_type> ()>
wrap_get_force(std::shared_ptr<particle_type> self)
{
    return [=]() -> std::vector<typename particle_type::force_type> {
        std::vector<typename particle_type::force_type> output;
        {
            output.reserve(self->force()->size());
        }
        get_force(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::en_pot_type> ()>
wrap_get_potential_energy(std::shared_ptr<particle_type> self)
{
    return [=]() -> std::vector<typename particle_type::en_pot_type> {
        std::vector<typename particle_type::en_pot_type> output;
        {
            output.reserve(self->potential_energy()->size());
        }
        get_potential_energy(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::stress_pot_type> ()>
wrap_get_stress_pot(std::shared_ptr<particle_type> self)
{
    return [=]() -> std::vector<typename particle_type::stress_pot_type> {
        std::vector<typename particle_type::stress_pot_type> output;
        {
            output.reserve(self->stress_pot()->size());
        }
        get_stress_pot(*self, std::back_inserter(output));
        return std::move(output);
    };
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
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
                    .def("get_position", &wrap_get_position<particle>)
                    .def("set_position", &wrap_set_position<particle>)
                    .def("get_image", &wrap_get_image<particle>)
                    .def("set_image", &wrap_set_image<particle>)
                    .def("get_velocity", &wrap_get_velocity<particle>)
                    .def("set_velocity", &wrap_set_velocity<particle>)
                    .def("get_tag", &wrap_get_tag<particle>)
                    .def("set_tag", &wrap_set_tag<particle>)
                    .def("get_reverse_tag", &wrap_get_reverse_tag<particle>)
                    .def("set_reverse_tag", &wrap_set_reverse_tag<particle>)
                    .def("get_species", &wrap_get_species<particle>)
                    .def("set_species", &wrap_set_species<particle>)
                    .def("get_mass", &wrap_get_mass<particle>)
                    .def("set_mass", &wrap_set_mass<particle>)
                    .def("get_force", &wrap_get_force<particle>)
                    .def("get_potential_energy", &wrap_get_potential_energy<particle>)
                    .def("get_stress_pot", &wrap_get_stress_pot<particle>)
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
