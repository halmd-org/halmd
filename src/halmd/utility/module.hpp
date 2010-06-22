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

#ifndef HALMD_UTILITY_MODULE_HPP
#define HALMD_UTILITY_MODULE_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/utility/modules/fetch.hpp>
#include <halmd/utility/modules/graph.hpp>
#include <halmd/utility/modules/traits.hpp>
#include <halmd/utility/modules/parser.hpp>
#include <halmd/utility/modules/registry.hpp>

namespace halmd
{

// shared_ptr is used in every module in HALMD, therefore
// we import it into the global halmd namespace.
using boost::shared_ptr;

namespace modules
{

template <typename Dependant, typename Dependency>
struct depends
{
    typedef modules::graph Graph;
    typedef modules::registry Registry;
    typedef typename boost::property_map<Graph, tag::relation>::type RelationPropertyMap;
    typedef typename boost::property_traits<RelationPropertyMap>::value_type RelationValue;
    typedef boost::color_traits<RelationValue> Relation;

    static void required()
    {
        Registry::template edge<Dependant, Dependency>(Relation::required());
    }

    static void optional()
    {
        Registry::template edge<Dependant, Dependency>(Relation::optional());
    }
};

} // namespace modules

template <typename T>
class module
{
public:
    typedef modules::typed_parser<T, modules::factory> Parser;

private:
    static Parser dummy_;
};

template <typename T>
typename module<T>::Parser module<T>::dummy_ = typename module<T>::Parser();

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULE_HPP */
