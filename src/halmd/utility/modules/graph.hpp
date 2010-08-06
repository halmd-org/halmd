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

#ifndef HALMD_UTILITY_MODULES_GRAPH_HPP
#define HALMD_UTILITY_MODULES_GRAPH_HPP

#define BOOST_NO_HASH  // circumvent deprecated header warning by disabling use of hash_set
#include <boost/graph/adjacency_list.hpp>
#undef BOOST_NO_HASH
#include <boost/graph/graph_traits.hpp>
#include <boost/shared_ptr.hpp>

#include <halmd/utility/modules/builder.hpp>
#include <halmd/utility/modules/property.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace modules
{

namespace tag
{

struct name         { typedef boost::vertex_property_tag kind; };
struct builder      { typedef boost::vertex_property_tag kind; };
struct selected     { typedef boost::vertex_property_tag kind; };
struct relation     { typedef boost::edge_property_tag kind; };

} // namespace tag

typedef boost::adjacency_list<
    boost::setS
  , boost::vecS
  , boost::bidirectionalS
  , boost::property<tag::name, std::string
      , boost::property<tag::builder, boost::shared_ptr<untyped_builder_base>
          , boost::property<tag::selected, boost::default_color_type
            >
        >
    >
  , boost::property<tag::relation, relation_type
    >
> graph;

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_GRAPH_HPP */
