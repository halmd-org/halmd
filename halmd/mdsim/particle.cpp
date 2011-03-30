/*
 * Copyright © 2008-2010  Peter Colberg
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
#include <boost/array.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <exception>
#include <numeric>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace boost::algorithm;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Construct microscopic system state.
 *
 * @param particles number of particles per type or species
 */
template <int dimension>
particle<dimension>::particle(vector<unsigned int> const& particles)
  : nbox(accumulate(particles.begin(), particles.end(), 0))
  , ntype(particles.size())
  , ntypes(particles)
{
    if (*min_element(this->ntypes.begin(), this->ntypes.end()) < 1) {
        throw logic_error("invalid number of particles");
    }

    vector<string> ntypes_(ntypes.size());
    transform(
        ntypes.begin()
      , ntypes.end()
      , ntypes_.begin()
      , lexical_cast<string, unsigned int>
    );

    LOG("number of particles: " << nbox);
    LOG("number of particle types: " << ntype);
    LOG("number of particles per type: " << join(ntypes_, " "));
}

template <int dimension>
void particle<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("particle_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("mdsim")
            [
                class_<particle, shared_ptr<particle> >(class_name.c_str())
                    .def_readonly("nbox", &particle::nbox)
                    .def_readonly("ntype", &particle::ntype)
                    .def_readonly("ntypes", &particle::ntypes)
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &particle<3>::luaopen
    ]
    [
        &particle<2>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class particle<3>;
template class particle<2>;

} // namespace mdsim

} // namespace halmd
