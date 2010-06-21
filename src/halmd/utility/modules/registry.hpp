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

#ifndef HALMD_UTILITY_MODULES_REGISTRY_HPP
#define HALMD_UTILITY_MODULES_REGISTRY_HPP

#include <halmd/utility/modules/graph.hpp>

namespace halmd
{
namespace modules
{

/**
 * The registry contains all singletons of the module mechanism,
 * i.e. the dependency graph and vertices for each module type.
 *
 * It is assembled before entrance into main() and shall not
 * be modified thereafter to guarantee thread safety.
 */
class registry
{
public:
    typedef modules::graph Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;

    //
    // Read this for the following use of static local variables:
    //
    // What's the "static initialization order fiasco"?
    // http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12
    //

    /**
     * This function constructs a global dependency graph
     * singleton that is filled before entrance into main().
     */
    static Graph& graph()
    {
        static Graph g;
        return g;
    }

    /**
     * Every module type is assigned a unique integer vertex in the
     * dependency graph; an simple type-to-vertex translation.
     */
    template <typename T>
    static Vertex const& vertex()
    {
        static Vertex v = _vertex<T>();
        return v;
    }

    template <typename Source, typename Target>
    static void edge(property::relation const& relation)
    {
        Graph& g = graph();
        Vertex u = vertex<Source>();
        Vertex v = vertex<Target>();
        add_edge(u, v, relation, g);
    }

private:
    /**
     * Create a new vertex and attach a type name property.
     */
    template <typename T>
    static Vertex _vertex()
    {
        Graph& g = graph();
        Vertex v = add_vertex(g);
        put(tag::name(), g, v, demangled_name<T>());
        put(tag::selected(), g, v, false);
        return v;
    }
};

}} // namespace halmd::modules

#endif /* ! HALMD_UTILITY_MODULES_REGISTRY_HPP */
