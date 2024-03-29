/*
 * Copyright © 2011  Peter Colberg
 * Copyright © 2012  Nicolas Höft
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

#include <algorithm>
#include <tuple>
#include <boost/bind/bind.hpp>

#include <halmd/algorithm/multi_range.hpp>
#include <halmd/mdsim/host/positions/excluded_volume.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace positions {

template <int dimension, typename float_type>
excluded_volume<dimension, float_type>::excluded_volume(
    std::shared_ptr<box_type const> box
  , float_type cell_length
  , std::shared_ptr<logger> logger
)
  : box_(box)
  , logger_(logger)
  , ncell_(vector_type(box_->length()) / cell_length)
  , cell_length_(element_div(vector_type(box_->length()), vector_type(ncell_)))
  , cell_(ncell_)
{
}

template <int dimension, typename float_type>
void excluded_volume<dimension, float_type>::exclude_sphere(
    vector_type const& centre
  , float_type diameter
)
{
    index_type lower, upper;
    std::tie(lower, upper) = this->sphere_extents(centre, diameter);
    upper += index_type(1); // index range is [lower, upper)
    multi_range_for_each(
        lower
      , upper
      , bind(&excluded_volume::exclude_sphere_from_cell, this, centre, diameter, boost::placeholders::_1)
    );
}

template <int dimension, typename float_type>
void excluded_volume<dimension, float_type>::exclude_spheres(
    position_sample_type const& position_sample
  , species_sample_type const& species_sample
  , std::vector<float_type> diameter
)
{
    typename position_sample_type::array_type const& position = position_sample.data();
    typename species_sample_type::array_type const& species = species_sample.data();

    for (size_t i = 0; i < position.size(); ++i) {
        vector_type r = position[i];
        unsigned int type = species[i];
        assert(type < diameter.size());
        this->exclude_sphere(r, diameter[type]);
    }
}

template <int dimension, typename float_type>
bool excluded_volume<dimension, float_type>::place_sphere(
    vector_type const& centre
  , float_type diameter
) const
{
    index_type lower, upper;
    std::tie(lower, upper) = this->sphere_extents(centre, diameter);
    upper += index_type(1); // index range is [lower, upper)
    index_type result = multi_range_find_if(
        lower
      , upper
      , bind(&excluded_volume::place_cell, this, centre, diameter, boost::placeholders::_1)
    );
    return std::equal(upper.begin(), upper.end(), result.begin());
}

template <int dimension, typename float_type>
typename excluded_volume<dimension, float_type>::index_pair_type
excluded_volume<dimension, float_type>::sphere_extents(
    vector_type const& centre
  , float_type diameter
) const
{
    vector_type lower = element_div(centre - vector_type(diameter / 2), cell_length_);
    vector_type upper = element_div(centre + vector_type(diameter / 2), cell_length_);
    for (unsigned int i = 0; i < dimension; ++i) {
        while (lower[i] < 0) {
            lower[i] += ncell_[i];
            upper[i] += ncell_[i];
        }
        assert(lower[i] >= 0);
        assert(upper[i] >= 0);
    }
    return std::make_pair(index_type(lower), index_type(upper));
}

template <int dimension, typename float_type>
void excluded_volume<dimension, float_type>::exclude_sphere_from_cell(
    vector_type const& centre
  , float_type diameter
  , index_type const& index
)
{
    // FIXME drop “spherical cow” approximation: a sphere is not a cube
    cell_(element_mod(index, ncell_)).push_back(std::make_pair(centre, diameter));
}

template <int dimension, typename float_type>
bool excluded_volume<dimension, float_type>::place_cell(
    vector_type const& centre
  , float_type diameter
  , index_type const& index
) const
{
    cell_type const& cell = cell_(element_mod(index, ncell_));
    for (size_t j = 0; j < cell.size(); ++j) {
        vector_type r = centre - cell[j].first;
        box_->reduce_periodic(r);
        float_type d = (diameter + cell[j].second) / 2;
        if (inner_prod(r, r) < d * d) {
            return true;
        }
    }
    return false;
}

template <int dimension, typename float_type>
void excluded_volume<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("positions")
            [
                class_<excluded_volume>()
                    .def("exclude_sphere", &excluded_volume::exclude_sphere)
                    .def("exclude_spheres", &excluded_volume::exclude_spheres)
                    .def("place_sphere", &excluded_volume::place_sphere)

              , def("excluded_volume", &std::make_shared<excluded_volume
                    , std::shared_ptr<box_type const>
                    , float_type
                    , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_positions_excluded_volume(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    excluded_volume<3, double>::luaopen(L);
    excluded_volume<2, double>::luaopen(L);
#else
    excluded_volume<3, float>::luaopen(L);
    excluded_volume<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class excluded_volume<3, double>;
template class excluded_volume<2, double>;
#else
template class excluded_volume<3, float>;
template class excluded_volume<2, float>;
#endif

} // namespace positions
} // namespace host
} // namespace mdsin
} // namespace halmd
