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

#ifndef HALMD_UTILITY_MODULES_PREDICATE_HPP
#define HALMD_UTILITY_MODULES_PREDICATE_HPP

#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/depth_first_search.hpp>

#include <halmd/utility/modules/graph.hpp>

namespace halmd
{
namespace modules { namespace predicate
{

/**
 * Given a relation property map and a relation property value,
 * this predicate includes only edges that match this value.
 */
template <typename PropertyMap>
struct relation
{
    PropertyMap map;
    property::relation pred;

    relation() {} // requirement of Iterator concept
    relation(PropertyMap const& map, property::relation const& pred)
      : map(map)
      , pred(pred)
    {}

    template <typename Edge>
    bool operator()(Edge const& e) const
    {
        return get(map, e) == pred;
    }
};

template <typename PropertyMap>
struct not_relation
{
    PropertyMap map;
    property::relation pred;

    not_relation() {} // requirement of Iterator concept
    not_relation(PropertyMap const& map, property::relation const& pred)
      : map(map)
      , pred(pred)
    {}

    template <typename Edge>
    bool operator()(Edge const& e) const
    {
        return get(map, e) != pred;
    }
};

/**
 * This predicate includes only vertices that are selected.
 */
template <typename PropertyMap>
struct selected
{
    PropertyMap map;

    selected() {} // requirement of Iterator concept
    explicit selected(PropertyMap const& map)
      : map(map)
    {}

    template <typename Vertex>
    bool operator()(Vertex const& v) const
    {
        return get(map, v); // boost::tribool == true
    }
};

/**
 * This predicate includes only vertices that are unsuitable.
 */
template <typename PropertyMap>
struct not_selected
{
    PropertyMap map;

    not_selected() {} // requirement of Iterator concept
    explicit not_selected(PropertyMap const& map)
      : map(map)
    {}

    template <typename Vertex>
    bool operator()(Vertex const& v) const
    {
        return !get(map, v); // boost::tribool == false
    }

    /**
     * DFS terminator
     */
    template <typename Vertex, typename Graph>
    bool operator()(Vertex const& v, Graph const&) const
    {
        return !get(map, v);
    }
};

template <typename PropertyMap>
struct not_not_selected
{
    PropertyMap map;

    not_not_selected() {} // requirement of Iterator concept
    explicit not_not_selected(PropertyMap const& map)
      : map(map)
    {}

    template <typename Vertex>
    bool operator()(Vertex const& v) const
    {
        return !bool(!get(map, v)); // Yikes!
    }
};

/**
 * This predicate includes only vertices without in-edges.
 */
template <typename Graph>
struct root
{
    Graph const* g; // use pointer to allow default constructor

    root() {} // requirement of Iterator concept
    explicit root(Graph const& g)
      : g(&g)
    {}

    template <typename Vertex>
    bool operator()(Vertex const& v) const
    {
        return !in_degree(v, *g);
    }
};

template <typename Graph>
struct selected_descendants
{
    Graph const* g; // use pointer to allow default constructor

    selected_descendants() {} // requirement of Iterator concept
    explicit selected_descendants(Graph const& g)
      : g(&g)
    {}

    template <typename Edge>
    bool operator()(Edge const& e) const
    {
        return get(tag::selected(), *g, target(e, *g));
    }
};

}} // namespace modules::predicate

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_PREDICATE_HPP */
