/*
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

#include <exception>

#include <halmd/algorithm/gpu/copy_if.hpp>
#include <halmd/mdsim/gpu/region.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <halmd/mdsim/geometries/cuboid.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * construct region module
 *
 * @param particle mdsim::gpu::particle instance
 * @param box mdsim::box instance
 */
template <int dimension, typename float_type, typename geometry_type>
region<dimension, float_type, geometry_type>::region(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<geometry_type> geometry
  , geometry_selection geometry_sel
  , std::shared_ptr<logger> logger
)
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  , geometry_(geometry)
  , geometry_selection_(geometry_sel)
{
    try {
        auto mask = make_cache_mutable(mask_);
        auto selection = make_cache_mutable(selection_);
        mask->reserve(particle_->dim().threads());
        mask->resize(particle_->nparticle());
        selection->reserve(particle_->dim().threads());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate global device memory");
        throw;
    }
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region<dimension, float_type, geometry_type>::array_type> const&
region<dimension, float_type, geometry_type>::selection()
{
    update_selection_();
    return selection_;
}

template <int dimension, typename float_type, typename geometry_type>
cache<typename region<dimension, float_type, geometry_type>::array_type> const&
region<dimension, float_type, geometry_type>::mask()
{
    update_mask_();
    return mask_;
}

/**
 * update the particle lists for the region
 */
template <int dimension, typename float_type, typename geometry_type>
void region<dimension, float_type, geometry_type>::update_mask_()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (position_cache != mask_cache_) {
        scoped_timer_type timer(runtime_.update_mask);

        auto mask = make_cache_mutable(mask_);
        position_array_type const& position = read_cache(particle_->position());
        auto const& kernel = region_wrapper<dimension, geometry_type>::kernel;
        // calculate "bin", ie. inside/outside the region
        cuda::memset(*mask, 0xFF);
        cuda::configure(particle_->dim().grid, particle_->dim().block);
        kernel.compute_mask(
            position.data()
          , particle_->nparticle()
          , mask->data()
          , *geometry_
          , geometry_selection_ == excluded ? halmd::mdsim::gpu::excluded : halmd::mdsim::gpu::included
          , static_cast<position_type>(box_->length())
        );
        mask_cache_ = position_cache;
    }
}

/**
 * update the particle lists for the region
 */
template <int dimension, typename float_type, typename geometry_type>
void region<dimension, float_type, geometry_type>::update_selection_()
{
    cache<position_array_type> const& position_cache = particle_->position();
    if(position_cache != selection_cache_) {
        scoped_timer_type timer(runtime_.update_selection);
        unsigned int nparticle = particle_->nparticle();
        auto const& position = read_cache(particle_->position());
        auto selection = make_cache_mutable(selection_);

        auto const& kernel = region_wrapper<dimension, geometry_type>::kernel;
        unsigned int size = kernel.copy_selection(
            position.data()
          , nparticle
          , selection->data()
          , *geometry_
          , geometry_selection_ == excluded ? halmd::mdsim::gpu::excluded : halmd::mdsim::gpu::included
        );
        selection->resize(size);

        selection_cache_ = position_cache;
    }
}

template <int dimension, typename float_type, typename geometry_type>
void region<dimension, float_type, geometry_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<region, region_base>()
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("update_mask", &runtime::update_mask)
                        .def_readonly("update_selection", &runtime::update_selection)
                ]
                .def_readonly("runtime", &region::runtime_)
          , def("region", &std::make_shared<region
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
                  , std::shared_ptr<geometry_type>
                  , geometry_selection
                  , std::shared_ptr<logger>
              >)
        ]
    ];
}

void region_base::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<region_base>()
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_region(lua_State* L)
{
    region_base::luaopen(L);
    region<3, float, mdsim::geometries::cuboid<3, float>>::luaopen(L);
    region<2, float, mdsim::geometries::cuboid<2, float>>::luaopen(L);
    return 0;
}

} // namespace gpu
} // namepsace mdsim
} // namespace halmd
