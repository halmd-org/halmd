/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/mdsim/host/neighbours/from_particle.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace neighbours {

/**
 * construct neighbour list module
 *
 * @param particle mdsim::host::particle instance
 * @param box mdsim::box instance
 * @param cutoff force cutoff radius
 * @param skin neighbour list skin
 */
template <int dimension, typename float_type>
from_particle<dimension, float_type>::from_particle(
    std::pair<std::shared_ptr<particle_type const>, std::shared_ptr<particle_type const>> particle
  , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>> displacement
  , std::shared_ptr<box_type const> box
  , matrix_type const& r_cut
  , double skin
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle1_(particle.first)
  , particle2_(particle.second)
  , displacement1_(displacement.first)
  , displacement2_(displacement.second)
  , box_(box)
  , logger_(logger)
  // allocate parameters
  , neighbour_(particle1_->nparticle())
  , r_skin_(skin)
  , rr_cut_skin_(particle1_->nspecies(), particle2_->nspecies())
{
    matrix_type r_cut_skin(r_cut.size1(), r_cut.size2());
    typename matrix_type::value_type r_cut_max = 0;
    for (size_t i = 0; i < r_cut.size1(); ++i) {
        for (size_t j = 0; j < r_cut.size2(); ++j) {
            r_cut_skin(i, j) = r_cut(i, j) + r_skin_;
            rr_cut_skin_(i, j) = std::pow(r_cut_skin(i, j), 2);
            r_cut_max = std::max(r_cut_skin(i, j), r_cut_max);
        }
    }

    LOG("neighbour list skin: " << r_skin_);
}

template <int dimension, typename float_type>
cache<std::vector<typename from_particle<dimension, float_type>::neighbour_list>> const&
from_particle<dimension, float_type>::lists()
{
    cache<reverse_tag_array_type> const& reverse_tag_cache1 = particle1_->reverse_tag();
    cache<reverse_tag_array_type> const& reverse_tag_cache2 = particle2_->reverse_tag();

    auto current_cache = std::tie(reverse_tag_cache1, reverse_tag_cache2);

    if (neighbour_cache_ != current_cache || displacement1_->compute() > r_skin_ / 2
        || displacement2_->compute() > r_skin_ / 2) {
        on_prepend_update_();
        update();
        displacement1_->zero();
        displacement2_->zero();
        neighbour_cache_ = current_cache;
        on_append_update_();
    }
    return neighbour_;
}

/**
 * Update neighbour lists
 */
template <int dimension, typename float_type>
void from_particle<dimension, float_type>::update()
{
    auto neighbour = make_cache_mutable(neighbour_);

    position_array_type const& position1 = read_cache(particle1_->position());
    position_array_type const& position2 = read_cache(particle2_->position());
    species_array_type const& species1 = read_cache(particle1_->species());
    species_array_type const& species2 = read_cache(particle2_->species());
    size_type nparticle1 = particle1_->nparticle();
    size_type nparticle2 = particle2_->nparticle();

    LOG_TRACE("update neighbour lists");

    scoped_timer_type timer(runtime_.update);

    // whether Newton's third law applies
    bool const reactio = (particle1_ == particle2_);

    for (size_type i = 0; i < nparticle1; ++i) {
        // load first particle
        vector_type r1 = position1[i];
        species_type type1 = species1[i];

        // clear particle's neighbour list
        (*neighbour)[i].clear();

        for (size_type j = reactio ? (i + 1) : 0; j < nparticle2; ++j) {
            // load second particle
            vector_type r2 = position2[j];
            species_type type2 = species2[j];

            // particle distance vector
            vector_type r = r1 - r2;
            box_->reduce_periodic(r);
            // squared particle distance
            float_type rr = inner_prod(r, r);

            // enforce cutoff radius with neighbour list skin
            if (rr >= rr_cut_skin_(type1, type2)) {
                continue;
            }

            // add particle to neighbour list
            (*neighbour)[i].push_back(j);
        }
    }
}

template <int dimension, typename float_type>
void from_particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("neighbours")
            [
                class_<from_particle, _Base>()
                    .property("r_skin", &from_particle::r_skin)
                    .def("on_prepend_update", &from_particle::on_prepend_update)
                    .def("on_append_update", &from_particle::on_append_update)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("update", &runtime::update)
                    ]
                    .def_readonly("runtime", &from_particle::runtime_)
              , def("from_particle", &std::make_shared<from_particle
                        , std::pair<std::shared_ptr<particle_type const>, std::shared_ptr<particle_type const>>
                        , std::pair<std::shared_ptr<displacement_type>, std::shared_ptr<displacement_type>>
                        , std::shared_ptr<box_type const>
                        , matrix_type const&
                        , double
                        , std::shared_ptr<logger>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_neighbours_from_particle(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    from_particle<3, double>::luaopen(L);
    from_particle<2, double>::luaopen(L);
#else
    from_particle<3, float>::luaopen(L);
    from_particle<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class from_particle<3, double>;
template class from_particle<2, double>;
#else
template class from_particle<3, float>;
template class from_particle<2, float>;
#endif

} // namespace neighbours
} // namespace host
} // namespace mdsim
} // namespace halmd
