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

#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/shared_ptr.hpp>

#include <halmd/utility/modules/builder.hpp>
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

namespace property
{

typedef std::string name;
typedef boost::shared_ptr<untyped_builder_base> builder;
typedef boost::default_color_type selected;
enum relation { is_required, is_optional, is_implicit, is_base_of };

} // namespace property

typedef boost::adjacency_list<
    boost::setS
  , boost::vecS
  , boost::bidirectionalS
  , boost::property<tag::name, property::name
      , boost::property<tag::builder, property::builder
          , boost::property<tag::selected, property::selected
            >
        >
    >
  , boost::property<tag::relation, property::relation
    >
> graph;

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_GRAPH_HPP */
