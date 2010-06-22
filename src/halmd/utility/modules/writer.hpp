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

#ifndef HALMD_UTILITY_MODULES_WRITER_HPP
#define HALMD_UTILITY_MODULES_WRITER_HPP

#include <boost/graph/graphviz.hpp>
#include <fstream>

#include <halmd/utility/modules/graph.hpp>

namespace halmd
{
namespace modules
{

/**
 * Format vertex properties for graphviz output.
 */
template <typename Graph>
struct vertex_property_writer
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::property_map<Graph, tag::name>::const_type NamePropertyMap;
    typedef typename boost::property_map<Graph, tag::builder>::const_type BuilderPropertyMap;
    typedef typename boost::property_map<Graph, tag::selected>::const_type SelectedPropertyMap;
    typedef typename boost::property_traits<SelectedPropertyMap>::value_type SelectedValue;
    typedef boost::color_traits<SelectedValue> Color;

    NamePropertyMap name;
    BuilderPropertyMap builder;
    SelectedPropertyMap selected;

    vertex_property_writer(Graph const& g)
      : name(get(tag::name(), g))
      , builder(get(tag::builder(), g))
      , selected(get(tag::selected(), g))
    {}

    void operator()(std::ostream& out, Vertex const& v) const
    {
        out << "[label=\"" << get(name, v) << "\"";
        if (get(builder, v)) {
            out << ",shape=\"box\"";
            if (get(selected, v) == Color::black()) {
                out << ",style=\"filled\",fillcolor=\"lightblue\"";
            }
            else if (get(selected, v) == Color::gray()) {
                out << ",style=\"filled\",fillcolor=\"lightgrey\"";
            }
            // else boost::indeterminate
        }
        else {
            out << ",shape=\"ellipse\"";
            if (get(selected, v) == Color::black()) {
                out << ",style=\"filled\",fillcolor=\"mistyrose\"";
            }
            else if (get(selected, v) == Color::gray()) {
                out << ",style=\"filled\",fillcolor=\"lightgrey\"";
            }
            // else boost::indeterminate
        }
        out << "]";
    }
};

template <typename Graph>
vertex_property_writer<Graph> make_vertex_property_writer(Graph const& g)
{
    return vertex_property_writer<Graph>(g);
}

/**
 * Format edge properties for graphviz output.
 */
template <typename Graph>
struct edge_property_writer
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::property_map<modules::graph, modules::tag::relation>::const_type RelationPropertyMap;
    typedef typename boost::property_traits<RelationPropertyMap>::value_type RelationValue;
    typedef boost::color_traits<RelationValue> Relation;

    RelationPropertyMap map;

    edge_property_writer(Graph const& g)
      : map(get(tag::relation(), g))
    {}

    void operator()(std::ostream& out, Edge const& e) const
    {
        RelationValue value = get(map, e);
        if (value == Relation::required())   out << "[style=\"solid\"]";
        if (value == Relation::optional())   out << "[style=\"dashed\"]";
        if (value == Relation::base())       out << "[arrowhead=\"ediamond\"]";
        if (value == Relation::implicit())   out << "[style=\"dotted\"]";
    }
};

template <typename Graph>
edge_property_writer<Graph> make_edge_property_writer(Graph const& g)
{
    return edge_property_writer<Graph>(g);
}

/**
 * Format graph properties for graphviz output.
 */
struct graph_property_writer
{
    void operator()(std::ostream& out) const
    {
        out << "graph [rankdir=\"LR\"]" << std::endl;
        out << "node [shape=\"box\"]" << std::endl;
        out << "edge [style=\"solid\"]" << std::endl;
    }
};

/**
 * Output a dependency graph of all modules in graphviz format.
 */
template <typename Graph>
void write_graphviz(std::string const& filename, Graph const& g)
{
    std::ofstream f;
    f.exceptions(std::ofstream::eofbit|std::ofstream::failbit|std::ofstream::badbit);
    f.open(filename.c_str());
    write_graphviz(
        f
      , g
      , make_vertex_property_writer(g)
      , make_edge_property_writer(g)
      , graph_property_writer()
    );
    f.close();
}

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_WRITER_HPP */
