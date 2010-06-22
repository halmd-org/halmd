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
        typedef typename boost::property_traits<RelationMap>::value_type RelationValue;
        typedef boost::color_traits<RelationValue> Relation;
        typedef predicate::relation<RelationMap> RelationPredicate;
        typedef boost::filtered_graph<Graph, RelationPredicate> RequiredGraph;
        typedef typename boost::graph_traits<RequiredGraph>::out_edge_iterator EdgeIterator;

        RelationPredicate ep(get(tag::relation(), g), Relation::required());
        RequiredGraph rg(g, ep);
        EdgeIterator ei, ei_end;
        for (boost::tie(ei, ei_end) = out_edges(source(e, g), rg); ei != ei_end; ++ei) {
            add_edge(target(e, g), target(*ei, g), Relation::implicit(), g);
        }
    }
};

template <typename Graph>
struct resolver
  : public boost::default_dfs_visitor
{
    typedef typename boost::property_map<Graph, tag::builder>::type BuilderPropertyMap;
    typedef typename boost::property_traits<BuilderPropertyMap>::value_type Builder;
    typedef typename boost::property_map<Graph, tag::selected>::type PropertyMap;
    typedef typename boost::property_traits<PropertyMap>::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;

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
        Builder builder = get(tag::builder(), g, v);
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
        put(tag::selected(), g, v, Color::gray());
    }

    template <typename Vertex, typename FilteredGraph>
    void finish_vertex(Vertex const& v, FilteredGraph const&)
    {
        typedef typename boost::property_map<Graph, tag::relation>::type RelationMap;
        typedef typename boost::property_traits<RelationMap>::value_type RelationValue;
        typedef boost::color_traits<RelationValue> Relation;
        typedef typename boost::property_map<Graph, tag::selected>::type SelectedMap;
        typedef predicate::relation<RelationMap> RelationPredicate;
        typedef predicate::selected<SelectedMap> SelectedPredicate;
        typedef predicate::not_selected<SelectedMap> NotSelectedPredicate;

        LOG_DEBUG("finish module " << get(tag::name(), g, v));
        if (get(tag::selected(), g, v) == Color::white()) {
            return;
        }
        RelationPredicate rp(get(tag::relation(), g), Relation::required());
        SelectedPredicate sp(get(tag::selected(), g), Color::white());
        if (out_degree(v, make_filtered_graph(g, rp, sp))) {
            LOG_DEBUG("✘ " << "missing required dependency");
            put(tag::selected(), g, v, Color::white());
            return;
        }
        if (!get(tag::builder(), g, v)) {
            RelationPredicate bp(get(tag::relation(), g), Relation::base());
            NotSelectedPredicate np(get(tag::selected(), g), Color::white());
            if (!out_degree(v, make_filtered_graph(g, bp, np))) {
                LOG_DEBUG("✘ " << "missing required module");
                put(tag::selected(), g, v, Color::white());
                return;
            }
        }
    }
};

template <typename Graph>
struct picker
  : public boost::default_dfs_visitor
{
    typedef typename boost::property_map<Graph, tag::selected>::type PropertyMap;
    typedef typename boost::property_traits<PropertyMap>::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef typename boost::property_map<Graph, tag::relation>::type RelationPropertyMap;
    typedef typename boost::property_traits<RelationPropertyMap>::value_type RelationValue;
    typedef boost::color_traits<RelationValue> Relation;

    Graph& g;
    picker(Graph& g) : g(g) {}

    template <typename Vertex, typename FilteredGraph>
    void start_vertex(Vertex const& v, FilteredGraph const&)
    {
        put(tag::selected(), g, v, Color::black());
    }

    template <typename Edge, typename FilteredGraph>
    void examine_edge(Edge const& e, FilteredGraph const&)
    {
        if (get(tag::relation(), g, e) != Relation::base()) {
            put(tag::selected(), g, target(e, g), Color::black());
        }
    }
};

template <typename Graph, typename RequiredStack>
struct policy
  : public boost::default_dfs_visitor
{
    typedef typename boost::property_map<Graph, tag::relation>::type RelationMap;
    typedef typename boost::property_map<Graph, tag::selected>::type SelectedMap;
    typedef typename boost::property_traits<SelectedMap>::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
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
                    put(tag::selected(), g, *ai, Color::black());
                }
            }
            put(tag::selected(), g, v, Color::gray());
        }
        stack.push_back(false);
    }

    template <typename Vertex, typename BaseFilteredGraph>
    void finish_vertex(Vertex const& v, BaseFilteredGraph const& bg)
    {
        if (get(tag::builder(), g, v)) {
            if (stack.back()) { // base module overriden by derived module
                LOG_DEBUG("✘ " << get(tag::name(), g, v));
                put(tag::selected(), g, v, Color::gray());
            }
            else if (get(tag::selected(), g, v) == Color::black()) {
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
        typedef typename boost::property_map<Graph, tag::selected>::type SelectedMap;
        typedef typename boost::property_traits<SelectedMap>::value_type ColorValue;
        typedef boost::color_traits<ColorValue> Color;
        typedef typename BuilderStack::iterator StackIterator;
        typedef typename BuilderMap::value_type::value_type Builder;

        stack.push_back(&map[v]);
        if (get(tag::selected(), g, v) == Color::black()) {
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
