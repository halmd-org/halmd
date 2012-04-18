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

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <exception>
#include <numeric>

#include <halmd/algorithm/host/permute.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {

/**
 * Allocate microscopic system state.
 *
 * @param particles number of particles per type or species
 */
template <int dimension, typename float_type>
particle<dimension, float_type>::particle(
    vector<unsigned int> const& particles
  , vector<double> const& mass
)
  : _Base(particles, mass)
  // allocate particle storage
  , r(nbox)
  , image(nbox)
  , v(nbox)
  , tag(nbox)
  , reverse_tag(nbox)
  , type(nbox)
  , force_(nbox)
  , en_pot_(nbox)
  , stress_pot_(nbox)
  , hypervirial_(nbox)
  // disable auxiliary variables by default
  , aux_flag_(false)
  , aux_valid_(false)
{
}

/**
 * set particle tags and types
 */
template <int dimension, typename float_type>
void particle<dimension, float_type>::set()
{
    // assign particle types
    for (size_t i = 0, j = 0; j < ntype; i += ntypes[j], ++j) {
        fill_n(type.begin() + i, ntypes[j], j);
    }
    // assign particle tags
    copy(counting_iterator<size_t>(0), counting_iterator<size_t>(nbox), tag.begin());
    // initially, the mapping tag → index is a 1:1 mapping
    copy(tag.begin(), tag.end(), reverse_tag.begin());
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

    algorithm::host::permute(r.begin(), r.end(), index.begin());
    algorithm::host::permute(image.begin(), image.end(), index.begin());
    algorithm::host::permute(v.begin(), v.end(), index.begin());
    // no permutation of forces
    algorithm::host::permute(tag.begin(), tag.end(), index.begin());
    algorithm::host::permute(type.begin(), type.end(), index.begin());

    // update reverse tags
    for (unsigned int i = 0; i < tag.size(); ++i) {
        assert(tag[i] < reverse_tag.size());
        reverse_tag[tag[i]] = i;
    }
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
}

template <typename particle_type>
static typename signal<void ()>::slot_function_type
wrap_aux_enable(shared_ptr<particle_type> self)
{
    return bind(&particle_type::aux_enable, self);
}

template <typename particle_type>
static typename signal<void ()>::slot_function_type
wrap_prepare(shared_ptr<particle_type> self)
{
    return bind(&particle_type::prepare, self);
}

template <typename particle_type>
struct wrap_particle
  : particle_type
  , luabind::wrap_base
{
    wrap_particle(
        vector<unsigned int> const& particles
      , vector<double> const& mass
    ) : particle_type(particles, mass) {}
};

template <int dimension, typename float_type>
void particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("particle_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<particle, shared_ptr<_Base>, _Base, wrap_particle<particle> >(class_name.c_str())
                    .def(constructor<
                        vector<unsigned int> const&
                      , vector<double> const&
                    >())
                    .property("nparticle", &particle::nparticle)
                    .property("nspecies", &particle::nspecies)
                    .property("dimension", &wrap_dimension<dimension, float_type>)
                    .property("aux_enable", &wrap_aux_enable<particle>)
                    .property("prepare", &wrap_prepare<particle>)
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

} // namespace mdsim
} // namespace host
} // namespace halmd
