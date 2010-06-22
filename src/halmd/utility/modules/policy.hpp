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

#ifndef HALMD_UTILITY_MODULES_POLICY_HPP
#define HALMD_UTILITY_MODULES_POLICY_HPP

#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/depth_first_search.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/modules/predicate.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace modules
{

class policy
{
public:
    typedef modules::registry Registry;
    typedef Registry::Graph Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef boost::default_color_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef std::vector<ColorValue> ColorMap;
    typedef std::vector<bool> RequiredStack;

    explicit policy(Graph const& g)
      : graph_(g)
    {
        typedef boost::property_map<Graph, tag::relation>::type RelationMap;
        typedef boost::property_map<Graph, tag::selected>::type SelectedMap;
        typedef predicate::relation<RelationMap> RelationPredicate;
        typedef predicate::selected<SelectedMap> SelectedPredicate;
        typedef predicate::not_selected<SelectedMap> NotSelectedPredicate;
        typedef boost::filtered_graph<Graph, RelationPredicate, SelectedPredicate> FilteredGraph;
        typedef boost::graph_traits<FilteredGraph>::edge_iterator EdgeIterator;
        typedef boost::graph_traits<FilteredGraph>::vertex_iterator VertexIterator;
        typedef predicate::root<FilteredGraph> RootPredicate;
        typedef boost::filtered_graph<FilteredGraph, boost::keep_all, RootPredicate> RootGraph;
        typedef boost::graph_traits<RootGraph>::vertex_iterator RootVertexIterator;
        typedef property::builder Builder;

        LOG_DEBUG("apply module policy");
        ColorMap color(num_vertices(graph_), Color::white());
        RequiredStack stack;
        RelationPredicate bp(get(tag::relation(), graph_), property::is_base_of);
        SelectedPredicate sp(get(tag::selected(), graph_), Color::black());
        NotSelectedPredicate np(get(tag::selected(), graph_), Color::white());
        FilteredGraph bg(graph_, bp, sp);
        RootGraph og(bg, boost::keep_all(), RootPredicate(bg));
        RootVertexIterator ri, ri_end;
        for (boost::tie(ri, ri_end) = vertices(og); ri != ri_end; ++ri) {
            depth_first_visit(
                make_filtered_graph(graph_, bp, np)
              , *ri // base class at bottom of class hierarchy
              , visitor::policy<Graph, RequiredStack>(graph_, stack)
              , &color.front()
            );
        }
    }

    Graph const& graph() const
    {
        return graph_;
    }

private:
    Graph graph_;
};

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_POLICY_HPP */
