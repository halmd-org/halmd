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

#ifndef HALMD_UTILITY_MODULES_RESOLVER_HPP
#define HALMD_UTILITY_MODULES_RESOLVER_HPP

#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/depth_first_search.hpp>

#include <halmd/utility/modules/predicate.hpp>
#include <halmd/utility/modules/registry.hpp>
#include <halmd/utility/modules/visitor.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace modules
{

class resolver
{
public:
    typedef modules::registry Registry;
    typedef Registry::Graph Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef boost::default_color_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef std::vector<ColorValue> ColorMap;

    explicit resolver(Graph const& g)
      : graph_(g)
      , color_(num_vertices(graph_), Color::white())
    {
        typedef boost::property_map<Graph, tag::relation>::type RelationMap;
        typedef predicate::relation<RelationMap> RelationPredicate;
        typedef boost::filtered_graph<Graph, RelationPredicate> FilteredGraph;
        typedef predicate::root<FilteredGraph> RootPredicate;
        typedef boost::filtered_graph<FilteredGraph, boost::keep_all, RootPredicate> RootGraph;
        typedef boost::graph_traits<RootGraph>::vertex_iterator VertexIterator;

        LOG_DEBUG("construct module resolver");
        RelationPredicate ep(get(tag::relation(), graph_), property::is_base_of);
        FilteredGraph fg(graph_, ep);
        RootGraph rg(fg, boost::keep_all(), RootPredicate(fg));
        VertexIterator vi, vi_end;
        ColorMap color(num_vertices(graph_), Color::white());
        for (boost::tie(vi, vi_end) = vertices(rg); vi != vi_end; ++vi) {
            depth_first_visit(
                fg
              , *vi // base class at bottom of class hierarchy
              , visitor::forwarder<Graph>(graph_)
              , &color.front()
            );
        }
    }

    template <typename T>
    void resolve(po::options const& vm, po::unparsed_options& unparsed)
    {
        typedef boost::property_map<Graph, tag::selected>::type SelectedMap;
        typedef predicate::not_selected<SelectedMap> NotSelectedPredicate;
        typedef predicate::not_not_selected<SelectedMap> NotNotSelectedPredicate;

        Vertex v = Registry::template vertex<T>();
        LOG_DEBUG("resolve module " << get(tag::name(), graph_, v));
        depth_first_visit(
            graph_
          , v
          , visitor::resolver<Graph>(graph_, vm, unparsed)
          , &color_.front()
          , NotSelectedPredicate(get(tag::selected(), graph_)) // terminate search
        );
        if (!get(tag::selected(), graph_, v)) {
            throw std::logic_error("failed to resolve module " + get(tag::name(), graph_, v));
        }
        NotNotSelectedPredicate nnp(get(tag::selected(), graph_));
        ColorMap color(num_vertices(graph_), Color::white());
        depth_first_visit(
            make_filtered_graph(graph_, boost::keep_all(), nnp)
          , v
          , visitor::picker<Graph>(graph_)
          , &color.front()
          , NotSelectedPredicate(get(tag::selected(), graph_)) // terminate search
        );
    }

    Graph const& graph() const
    {
        return graph_;
    }

private:
    Graph graph_;
    ColorMap color_;
};

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_RESOLVER_HPP */
