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

#include <halmd/config.hpp>

#include <algorithm>
#include <boost/function_output_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <exception>
#include <iterator> // std::back_inserter
#include <luabind/luabind.hpp>
#include <luabind/out_value_policy.hpp>
#include <numeric>

#include <halmd/algorithm/host/permute.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/particle_group.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
particle<dimension, float_type>::particle(size_t nparticle, unsigned int nspecies)
  // allocate particle storage
  : nspecies_(std::max(nspecies, 1u))
  , position_(nparticle, 0)
  , image_(nparticle, 0)
  , velocity_(nparticle, 0)
  , tag_(counting_iterator<tag_type>(0), counting_iterator<tag_type>(nparticle))
  , reverse_tag_(tag_)
  , species_(nparticle, 0)
  , mass_(nparticle, 1)
  , force_(nparticle, 0)
  , en_pot_(nparticle, 0)
  , stress_pot_(nparticle, 0)
  , hypervirial_(nparticle, 0)
  // disable auxiliary variables by default
  , aux_flag_(false)
  , aux_valid_(false)
{
    LOG("number of particles: " << tag_.size());
    LOG("number of particle placeholders: " << tag_.capacity());
    LOG("number of particle species: " << nspecies_);
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::set_mass(float_type mass)
{
    fill(mass_.begin(), mass_.end(), mass);
}

/**
 * set particle tags and types
 */
template <int dimension, typename float_type>
void particle<dimension, float_type>::set()
{
    // initialise particle types to zero
    fill(species_.begin(), species_.end(), 0);
    // assign particle tags
    copy(counting_iterator<size_t>(0), counting_iterator<size_t>(tag_.size()), tag_.begin());
    // initially, the mapping tag → index is a 1:1 mapping
    copy(tag_.begin(), tag_.end(), reverse_tag_.begin());
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::aux_enable()
{
    LOG_TRACE("enable computation of auxiliary variables");
    aux_flag_ = true;
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::prepare()
{
    LOG_TRACE("zero forces");
    fill(force_.begin(), force_.end(), 0);

    // indicate whether auxiliary variables are computed this step
    aux_valid_ = aux_flag_;

    if (aux_flag_) {
        LOG_TRACE("zero auxiliary variables");
        fill(en_pot_.begin(), en_pot_.end(), 0);
        fill(stress_pot_.begin(), stress_pot_.end(), 0);
        fill(hypervirial_.begin(), hypervirial_.end(), 0);
        aux_flag_ = false;
    }
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

    algorithm::host::permute(position_.begin(), position_.end(), index.begin());
    algorithm::host::permute(image_.begin(), image_.end(), index.begin());
    algorithm::host::permute(velocity_.begin(), velocity_.end(), index.begin());
    // no permutation of forces
    algorithm::host::permute(tag_.begin(), tag_.end(), index.begin());
    algorithm::host::permute(species_.begin(), species_.end(), index.begin());
    algorithm::host::permute(mass_.begin(), mass_.end(), index.begin());

    // update reverse tags
    for (unsigned int i = 0; i < tag_.size(); ++i) {
        assert(tag_[i] < reverse_tag_.size());
        reverse_tag_[tag_[i]] = i;
    }
}

template <typename particle_type>
static void
wrap_get_position(particle_type const& particle, vector<typename particle_type::position_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_position(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_position(particle_type& particle, vector<typename particle_type::position_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_position(input.begin(), input.end());
}

template <typename particle_type>
static void
wrap_get_image(particle_type const& particle, vector<typename particle_type::image_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_image(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_image(particle_type& particle, vector<typename particle_type::image_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_image(input.begin(), input.end());
}

template <typename particle_type>
static void
wrap_get_velocity(particle_type const& particle, vector<typename particle_type::velocity_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_velocity(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_velocity(particle_type& particle, vector<typename particle_type::velocity_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_velocity(input.begin(), input.end());
}

template <typename particle_type>
static void
wrap_get_tag(particle_type const& particle, vector<typename particle_type::tag_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    auto convert_0_to_1 = [&](typename particle_type::tag_type t) {
        output.push_back(t + 1);
    };
    particle.get_tag(
        boost::make_function_output_iterator(convert_0_to_1)
    );
}

template <typename particle_type>
static void
wrap_set_tag(particle_type& particle, vector<typename particle_type::tag_type> const& input)
{
    typedef typename particle_type::tag_type tag_type;
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    tag_type nparticle = particle.nparticle();
    auto convert_1_to_0 = [&](tag_type t) -> tag_type {
        if (t < 1 || t > nparticle) {
            throw std::invalid_argument("invalid particle tag");
        }
        return t - 1;
    };
    particle.set_tag(
        boost::make_transform_iterator(input.begin(), convert_1_to_0)
      , boost::make_transform_iterator(input.end(), convert_1_to_0)
    );
}

template <typename particle_type>
static void
wrap_get_reverse_tag(particle_type const& particle, vector<typename particle_type::reverse_tag_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    auto convert_0_to_1 = [&](typename particle_type::reverse_tag_type i) {
        output.push_back(i + 1);
    };
    particle.get_reverse_tag(
        boost::make_function_output_iterator(convert_0_to_1)
    );
}

template <typename particle_type>
static void
wrap_set_reverse_tag(particle_type& particle, vector<typename particle_type::reverse_tag_type> const& input)
{
    typedef typename particle_type::reverse_tag_type reverse_tag_type;
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    reverse_tag_type nparticle = particle.nparticle();
    auto convert_1_to_0 = [&](reverse_tag_type i) -> reverse_tag_type {
        if (i < 1 || i > nparticle) {
            throw std::invalid_argument("invalid particle reverse tag");
        }
        return i - 1;
    };
    particle.set_reverse_tag(
        boost::make_transform_iterator(input.begin(), convert_1_to_0)
      , boost::make_transform_iterator(input.end(), convert_1_to_0)
    );
}

template <typename particle_type>
static void
wrap_get_species(particle_type const& particle, vector<typename particle_type::species_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    auto convert_0_to_1 = [&](typename particle_type::species_type s) {
        output.push_back(s + 1);
    };
    particle.get_species(
        boost::make_function_output_iterator(convert_0_to_1)
    );
}

template <typename particle_type>
static void
wrap_set_species(particle_type& particle, vector<typename particle_type::species_type> const& input)
{
    typedef typename particle_type::species_type species_type;
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    species_type nspecies = particle.nspecies();
    auto convert_1_to_0 = [&](species_type s) -> species_type {
        if (s < 1 || s > nspecies) {
            throw std::invalid_argument("invalid particle species");
        }
        return s - 1;
    };
    particle.set_species(
        boost::make_transform_iterator(input.begin(), convert_1_to_0)
      , boost::make_transform_iterator(input.end(), convert_1_to_0)
    );
}

template <typename particle_type>
static void
wrap_get_mass(particle_type const& particle, vector<typename particle_type::mass_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_mass(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_mass(particle_type& particle, vector<typename particle_type::mass_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_mass(input.begin(), input.end());
}

template <typename particle_type>
static void
wrap_get_force(particle_type const& particle, vector<typename particle_type::force_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_force(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_force(particle_type& particle, vector<typename particle_type::force_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_force(input.begin(), input.end());
}

template <typename particle_type>
static void
wrap_get_en_pot(particle_type const& particle, vector<typename particle_type::en_pot_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_en_pot(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_en_pot(particle_type& particle, vector<typename particle_type::en_pot_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_en_pot(input.begin(), input.end());
}

template <typename particle_type>
static void
wrap_get_stress_pot(particle_type const& particle, vector<typename particle_type::stress_pot_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_stress_pot(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_stress_pot(particle_type& particle, vector<typename particle_type::stress_pot_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_stress_pot(input.begin(), input.end());
}

template <typename particle_type>
static void
wrap_get_hypervirial(particle_type const& particle, vector<typename particle_type::hypervirial_type>& output)
{
    output.clear();
    output.reserve(particle.nparticle());
    particle.get_hypervirial(back_inserter(output));
}

template <typename particle_type>
static void
wrap_set_hypervirial(particle_type& particle, vector<typename particle_type::hypervirial_type> const& input)
{
    if (input.size() != particle.nparticle()) {
        throw invalid_argument("input array size not equal to number of particles");
    }
    particle.set_hypervirial(input.begin(), input.end());
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
}

template <typename particle_type>
static typename signal<void ()>::slot_function_type
wrap_aux_enable(boost::shared_ptr<particle_type> self)
{
    return bind(&particle_type::aux_enable, self);
}

template <typename particle_type>
static typename signal<void ()>::slot_function_type
wrap_prepare(boost::shared_ptr<particle_type> self)
{
    return bind(&particle_type::prepare, self);
}

template <typename particle_type>
typename signal<void ()>::slot_function_type
wrap_set(boost::shared_ptr<particle_type> particle)
{
    return bind(&particle_type::set, particle);
}

template <typename particle_type>
struct wrap_particle
  : particle_type
  , luabind::wrap_base
{
    wrap_particle(size_t nparticle, unsigned int nspecies) : particle_type(nparticle, nspecies) {}
};

template <int dimension, typename float_type>
void particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("particle_" + lexical_cast<string>(dimension));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<particle, boost::shared_ptr<particle>, wrap_particle<particle> >(class_name.c_str())
                    .def(constructor<size_t, unsigned int>())
                    .property("nparticle", &particle::nparticle)
                    .property("nspecies", &particle::nspecies)
                    .def("get_position", &wrap_get_position<particle>, pure_out_value(_2))
                    .def("set_position", &wrap_set_position<particle>)
                    .def("get_image", &wrap_get_image<particle>, pure_out_value(_2))
                    .def("set_image", &wrap_set_image<particle>)
                    .def("get_velocity", &wrap_get_velocity<particle>, pure_out_value(_2))
                    .def("set_velocity", &wrap_set_velocity<particle>)
                    .def("get_tag", &wrap_get_tag<particle>, pure_out_value(_2))
                    .def("set_tag", &wrap_set_tag<particle>)
                    .def("get_reverse_tag", &wrap_get_reverse_tag<particle>, pure_out_value(_2))
                    .def("set_reverse_tag", &wrap_set_reverse_tag<particle>)
                    .def("get_species", &wrap_get_species<particle>, pure_out_value(_2))
                    .def("set_species", &wrap_set_species<particle>)
                    .def("get_mass", &wrap_get_mass<particle>, pure_out_value(_2))
                    .def("set_mass", &wrap_set_mass<particle>)
                    .def("set_mass", static_cast<void (particle::*)(float_type)>(&particle::set_mass))
                    .def("get_force", &wrap_get_force<particle>, pure_out_value(_2))
                    .def("set_force", &wrap_set_force<particle>)
                    .def("get_en_pot", &wrap_get_en_pot<particle>, pure_out_value(_2))
                    .def("set_en_pot", &wrap_set_en_pot<particle>)
                    .def("get_stress_pot", &wrap_get_stress_pot<particle>, pure_out_value(_2))
                    .def("set_stress_pot", &wrap_set_stress_pot<particle>)
                    .def("get_hypervirial", &wrap_get_hypervirial<particle>, pure_out_value(_2))
                    .def("set_hypervirial", &wrap_set_hypervirial<particle>)
                    .property("dimension", &wrap_dimension<dimension, float_type>)
                    .property("aux_enable", &wrap_aux_enable<particle>)
                    .property("prepare", &wrap_prepare<particle>)
                    .property("set", &wrap_set<particle>)
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
    particle_group<particle<3, double> >::luaopen(L);
    particle_group<particle<2, double> >::luaopen(L);
    particle_group_from_range<particle<3, double> >::luaopen(L);
    particle_group_from_range<particle<2, double> >::luaopen(L);
#else
    particle<3, float>::luaopen(L);
    particle<2, float>::luaopen(L);
    particle_group<particle<3, float> >::luaopen(L);
    particle_group<particle<2, float> >::luaopen(L);
    particle_group_from_range<particle<3, float> >::luaopen(L);
    particle_group_from_range<particle<2, float> >::luaopen(L);
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

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class particle_group_from_range<host::particle<3, double> >;
template class particle_group_from_range<host::particle<2, double> >;
#else
template class particle_group_from_range<host::particle<3, float> >;
template class particle_group_from_range<host::particle<2, float> >;
#endif

} // namespace mdsim
} // namespace halmd
