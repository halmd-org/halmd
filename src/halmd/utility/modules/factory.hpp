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

#ifndef HALMD_UTILITY_MODULES_FACTORY_HPP
#define HALMD_UTILITY_MODULES_FACTORY_HPP

#include <boost/graph/filtered_graph.hpp>

#include <halmd/utility/modules/predicate.hpp>
#include <halmd/utility/modules/registry.hpp>
#include <halmd/utility/modules/visitor.hpp>

namespace halmd
{
namespace modules
{

class factory
{
public:
    typedef modules::registry Registry;
    typedef Registry::Graph Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef boost::default_color_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef std::vector<ColorValue> ColorMap;
    typedef property::builder Builder;
    typedef std::map<Vertex, std::vector<Builder> > BuilderMap;

    explicit factory(Graph const& g)
    {
        typedef boost::property_map<Graph, tag::relation>::const_type RelationMap;
        typedef boost::property_map<Graph, tag::selected>::const_type SelectedMap;
        typedef predicate::relation<RelationMap> RelationPredicate;
        typedef predicate::not_selected<SelectedMap> NotSelectedPredicate;
        typedef boost::filtered_graph<Graph, RelationPredicate, NotSelectedPredicate> FilteredGraph;
        typedef predicate::root<FilteredGraph> RootPredicate;
        typedef boost::filtered_graph<FilteredGraph, boost::keep_all, RootPredicate> RootGraph;
        typedef boost::graph_traits<RootGraph>::vertex_iterator VertexIterator;
        typedef std::vector<std::vector<Builder>*> BuilderStack;

        LOG_DEBUG("construct module factory");
        BuilderStack stack;
        ColorMap color(num_vertices(g), Color::white()); // manually color due to partial DFS
        RelationPredicate ep(get(tag::relation(), g), property::is_base_of);
        NotSelectedPredicate np(get(tag::selected(), g), Color::white());
        FilteredGraph fg(g, ep, np);
        RootGraph rg(fg, boost::keep_all(), RootPredicate(fg));
        VertexIterator vi, vi_end;
        for (boost::tie(vi, vi_end) = vertices(rg); vi != vi_end; ++vi) {
            depth_first_visit(
                fg
              , *vi // base class at bottom of class hierarchy
              , visitor::factory<BuilderMap, BuilderStack>(builder, stack)
              , &color.front()
            );
        }
    }

    BuilderMap builder;
};

}} // namespace halmd::modules

#endif /* ! HALMD_UTILITY_MODULES_FACTORY_HPP */
