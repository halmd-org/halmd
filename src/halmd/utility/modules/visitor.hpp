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

#ifndef HALMD_UTILITY_MODULES_VISITOR_HPP
#define HALMD_UTILITY_MODULES_VISITOR_HPP

#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/depth_first_search.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/modules/predicate.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace modules { namespace visitor
{

template <typename Graph>
struct forwarder
  : public boost::default_dfs_visitor
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

    forwarder(Graph& g) : g(g) {}
    Graph& g;

    template <typename AcyclicGraph>
    void examine_edge(Edge const& e, AcyclicGraph const&)
    {
        typedef typename boost::property_map<Graph, tag::relation>::type RelationMap;
        typedef predicate::relation<RelationMap> RelationPredicate;
        typedef boost::filtered_graph<Graph, RelationPredicate> RequiredGraph;
        typedef typename boost::graph_traits<RequiredGraph>::out_edge_iterator EdgeIterator;

        RelationPredicate ep(get(tag::relation(), g), property::is_required);
        RequiredGraph rg(g, ep);
        EdgeIterator ei, ei_end;
        for (boost::tie(ei, ei_end) = out_edges(source(e, g), rg); ei != ei_end; ++ei) {
            add_edge(target(e, g), target(*ei, g), property::is_implicit, g);
        }
    }
};

template <typename Graph>
struct resolver
  : public boost::default_dfs_visitor
{
    resolver(Graph& g, po::options const& vm, po::unparsed_options& unparsed)
      : g(g)
      , vm(vm)
      , unparsed(unparsed)
    {}

    Graph& g;
    po::options const& vm;
    po::unparsed_options& unparsed;

    template <typename Vertex, typename FilteredGraph>
    void discover_vertex(Vertex const& v, FilteredGraph const&)
    {
        LOG_DEBUG("discover module " << get(tag::name(), g, v));
        property::builder builder = get(tag::builder(), g, v);
        if (builder) {
            builder->vm = vm;
            po::options_description desc;
            builder->options(desc);
            try {
                po::parse_options(unparsed, desc, builder->vm);
            }
            catch (po::required_option const& e) {
                LOG_DEBUG("✘ " << e.what());
                return;
            }
            try {
                builder->select(builder->vm);
            }
            catch (module_error const& e) {
                LOG_DEBUG("✘ " << e.what());
                return;
            }
        }
        put(tag::selected(), g, v, boost::tribool(boost::indeterminate));
    }

    template <typename Vertex, typename FilteredGraph>
    void finish_vertex(Vertex const& v, FilteredGraph const&)
    {
        typedef typename boost::property_map<Graph, tag::relation>::type RelationMap;
        typedef typename boost::property_map<Graph, tag::selected>::type SelectedMap;
        typedef predicate::relation<RelationMap> RelationPredicate;
        typedef predicate::not_selected<SelectedMap> NotSelectedPredicate;
        typedef predicate::not_not_selected<SelectedMap> NotNotSelectedPredicate;

        LOG_DEBUG("finish module " << get(tag::name(), g, v));
        if (!get(tag::selected(), g, v)) {
            return;
        }
        RelationPredicate rp(get(tag::relation(), g), property::is_required);
        NotSelectedPredicate np(get(tag::selected(), g));
        if (out_degree(v, make_filtered_graph(g, rp, np))) {
            LOG_DEBUG("✘ " << "missing required dependency");
            put(tag::selected(), g, v, false);
            return;
        }
        if (!get(tag::builder(), g, v)) {
            RelationPredicate bp(get(tag::relation(), g), property::is_base_of);
            NotNotSelectedPredicate nnp(get(tag::selected(), g));
            if (!out_degree(v, make_filtered_graph(g, bp, nnp))) {
                LOG_DEBUG("✘ " << "missing required module");
                put(tag::selected(), g, v, false);
                return;
            }
        }
    }
};

template <typename Graph>
struct picker
  : public boost::default_dfs_visitor
{
    Graph& g;
    picker(Graph& g) : g(g) {}

    template <typename Vertex, typename FilteredGraph>
    void start_vertex(Vertex const& v, FilteredGraph const&)
    {
        put(tag::selected(), g, v, true);
    }

    template <typename Edge, typename FilteredGraph>
    void examine_edge(Edge const& e, FilteredGraph const&)
    {
        if (get(tag::relation(), g, e) != property::is_base_of) {
            put(tag::selected(), g, target(e, g), true);
        }
    }
};

template <typename Graph, typename RequiredStack>
struct policy
  : public boost::default_dfs_visitor
{
    typedef typename boost::property_map<Graph, tag::relation>::type RelationMap;
    typedef typename boost::property_map<Graph, tag::selected>::type SelectedMap;
    typedef predicate::not_relation<RelationMap> NotRelationPredicate;
    typedef predicate::selected<SelectedMap> SelectedPredicate;
    typedef boost::filtered_graph<Graph, NotRelationPredicate, SelectedPredicate> FilteredGraph;
    typedef predicate::selected_descendants<Graph> SelectedDescendantsPredicate;

    Graph& g;
    RequiredStack& stack;

    policy(Graph& g, RequiredStack& stack)
      : g(g)
      , stack(stack)
    {}

    template <typename Vertex, typename BaseFilteredGraph>
    void discover_vertex(Vertex const& v, BaseFilteredGraph const& bg)
    {
        typedef typename boost::graph_traits<BaseFilteredGraph>::adjacency_iterator AdjacencyIterator;

        if (!get(tag::builder(), g, v)) {
            if (!out_degree(v, make_filtered_graph(bg, SelectedDescendantsPredicate(g)))) {
                AdjacencyIterator ai, ai_end;
                for (boost::tie(ai, ai_end) = adjacent_vertices(v, bg); ai != ai_end; ++ai) {
                    put(tag::selected(), g, *ai, true);
                }
            }
            put(tag::selected(), g, v, boost::tribool(boost::indeterminate));
        }
        stack.push_back(false);
    }

    template <typename Vertex, typename BaseFilteredGraph>
    void finish_vertex(Vertex const& v, BaseFilteredGraph const& bg)
    {
        if (get(tag::builder(), g, v)) {
            if (stack.back()) { // base module overriden by derived module
                LOG_DEBUG("✘ " << get(tag::name(), g, v));
                put(tag::selected(), g, v, boost::tribool(boost::indeterminate));
            }
            else if (get(tag::selected(), g, v)) {
                std::fill(stack.begin(), stack.end(), true);
            }
        }
        stack.pop_back();
    }
};

template <typename BuilderMap, typename BuilderStack>
struct factory
  : public boost::default_dfs_visitor
{
    BuilderMap& map;
    BuilderStack& stack;

    factory(BuilderMap& map, BuilderStack& stack)
      : map(map)
      , stack(stack)
    {}

    template <typename Vertex, typename Graph>
    void discover_vertex(Vertex const& v, Graph const& g)
    {
        typedef typename BuilderStack::iterator StackIterator;
        typedef typename BuilderMap::value_type::second_type::value_type Builder;

        stack.push_back(&map[v]);
        if (get(tag::selected(), g, v)) {
            Builder builder = get(tag::builder(), g, v);
            if (builder) {
                LOG_DEBUG("✔ " << get(tag::name(), g, v));
                StackIterator si, si_end;
                for (boost::tie(si, si_end) = std::make_pair(stack.begin(), stack.end()); si != si_end; ++si) {
                    (*si)->push_back(builder);
                }
            }
        }
    }

    template <typename Vertex, typename Graph>
    void finish_vertex(Vertex const& v, Graph const& g)
    {
        stack.pop_back();
    }
};

}} // namespace modules::visitor

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_VISITOR_HPP */
