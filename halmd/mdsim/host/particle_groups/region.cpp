/* Copyright © 2019 Roya Ebrahimi Viand
 * Copyright © 2014-2015 Nicolas Höft
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

#include <halmd/mdsim/geometries/cuboid.hpp>
#include <halmd/mdsim/geometries/sphere.hpp>
#include <halmd/algorithm/host/radix_sort.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/region.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <algorithm>

namespace halmd {
namespace mdsim {
namespace host {
namespace particle_groups {

template <int dimension, typename float_type, typename geometry_type>
region<dimension, float_type, geometry_type>::region(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<geometry_type const> geometry
  , geometry_selection geometry_sel
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  , geometry_(geometry)
  , geometry_selection_(geometry_sel)
  , mask_(particle->nparticle())
{
}

/**
 * update the particle lists for the region
 */
template <int dimension, typename float_type, typename geometry_type>
void region<dimension, float_type, geometry_type>::update_()
{
    cache<position_array_type> const& position_cache = particle_->position();
    if (position_cache != mask_cache_) {
        LOG_TRACE("update region");

        scoped_timer_type timer(runtime_.update_mask);

        auto mask = make_cache_mutable(mask_);
        auto selection = make_cache_mutable(selection_);
        position_array_type const& position = read_cache(particle_->position());

        selection->clear();

        for (size_type i = 0; i < particle_->nparticle(); ++i) {
            vector_type r = position[i];
            // box_->reduce_periodic(r);  FIXME test whether this is really not needed
            bool in_geometry  = (*geometry_)(r);
            if (geometry_selection_ == excluded) {
                in_geometry = !in_geometry;
            }
            (*mask)[i] = in_geometry ? 1 : 0;
            if (in_geometry) {
                selection->push_back(i);
            }
        }
        mask_cache_ = position_cache;
    }
}
template <int dimension, typename float_type, typename geometry_type>
cache<typename region<dimension, float_type, geometry_type>::array_type> const&
region<dimension, float_type, geometry_type>::selection()
{
    update_();
    return selection_;
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region<dimension, float_type, geometry_type>::array_type> const&
region<dimension, float_type, geometry_type>::mask()
{
    update_();
    return mask_;
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region<dimension, float_type, geometry_type>::array_type> const&
region<dimension, float_type, geometry_type>::ordered()
{
    auto const& selection_cache = selection();
    if (selection_cache != ordered_cache_) {
        auto ordered = make_cache_mutable(ordered_);
        ordered->clear(); // avoid copying the elements upon resize()
        ordered->resize(selection_cache->size());
        std::copy(selection_cache->begin(), selection_cache->end(), ordered->begin());

        ordered_cache_ = selection_cache;
    }
    return ordered_;
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region<dimension, float_type, geometry_type>::array_type> const&
region<dimension, float_type, geometry_type>::unordered()
{
    auto const& selection_cache = selection();
    if (selection_cache != unordered_cache_) {
        auto unordered = make_cache_mutable(unordered_);
        LOG_TRACE("unordered sequence of particle indices");

        unordered->clear(); // avoid copying the elements upon resize()
        unordered->resize(selection_cache->size());
            std::copy(selection_cache->begin(), selection_cache->end(), unordered->begin());

            // TODO: is radix sort required here?
            radix_sort(
                unordered->begin()
              , unordered->end()
            );

            unordered_cache_ = selection_cache;
        }
        return unordered_;
    }

template <int dimension, typename float_type, typename geometry_type>
cache<typename region<dimension, float_type, geometry_type>::size_type> const&
region<dimension, float_type, geometry_type>::size()
{
    auto const& s = selection();
    if (s != size_cache_) {
        auto size = make_cache_mutable(size_);
        *size = selection()->size();
        size_cache_ = s;
    }
    return size_;
}

    template <typename particle_group_type, typename particle_type>
    static void
    wrap_to_particle(std::shared_ptr<particle_group_type> self, std::shared_ptr<particle_type> particle_src, std::shared_ptr<particle_type> particle_dst)
    {
        particle_group_to_particle(*particle_src, *self, *particle_dst);
    }

    template <int dimension, typename float_type, typename geometry_type>
    void region<dimension, float_type, geometry_type>::luaopen(lua_State* L)
    {
        using namespace luaponte;
        module(L, "libhalmd")
        [
            namespace_("mdsim")
            [
                namespace_("particle_groups")
                [
                    class_<region, particle_group>()
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("update_mask", &runtime::update_mask)
                                .def_readonly("update_selection", &runtime::update_selection)
                        ]
                        .def_readonly("runtime", &region::runtime_)
                        .def("to_particle", &wrap_to_particle<region<dimension, float_type, geometry_type>, particle_type>)

                  , def("region", &std::make_shared<region<dimension, float_type, geometry_type>
                      , std::shared_ptr<particle_type const>
                      , std::shared_ptr<box_type const>
                      , std::shared_ptr<geometry_type const>
                      , geometry_selection
                      , std::shared_ptr<logger>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle_groups_region(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    region<3, double, halmd::mdsim::geometries::cuboid<3, double>>::luaopen(L);
    region<2, double, halmd::mdsim::geometries::cuboid<2, double>>::luaopen(L);
    region<3, double, halmd::mdsim::geometries::sphere<3, double>>::luaopen(L);
    region<2, double, halmd::mdsim::geometries::sphere<2, double>>::luaopen(L);
#else
    region<3, float, halmd::mdsim::geometries::cuboid<3, float>>::luaopen(L);
    region<2, float, halmd::mdsim::geometries::cuboid<2, float>>::luaopen(L);
    region<3, float, halmd::mdsim::geometries::sphere<3, float>>::luaopen(L);
    region<2, float, halmd::mdsim::geometries::sphere<2, float>>::luaopen(L);
#endif
    return 0;
}

} // namespace particle_groups
} // namespace host
} // namespace mdsim
} // namespace halmd
