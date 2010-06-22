/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_UTILITY_MODULES_PROPERTY_HPP
#define HALMD_UTILITY_MODULES_PROPERTY_HPP

#include <halmd/utility/modules/graph.hpp>

namespace halmd
{
namespace modules
{

enum relation_type { required_relation, optional_relation, implicit_relation, base_relation };

} // namespace modules

} // namespace halmd

namespace boost
{

template <>
struct color_traits<halmd::modules::relation_type>
{
    typedef halmd::modules::relation_type RelationValue;
    static RelationValue required() { return halmd::modules::required_relation; }
    static RelationValue optional() { return halmd::modules::optional_relation; }
    static RelationValue implicit() { return halmd::modules::implicit_relation; }
    static RelationValue base() { return halmd::modules::base_relation; }
};

} // namespace boost

#endif /* ! HALMD_UTILITY_MODULES_PROPERTY_HPP */
