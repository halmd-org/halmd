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

#ifndef HALMD_UTILITY_MODULES_PARSER_HPP
#define HALMD_UTILITY_MODULES_PARSER_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/type_traits/is_object.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/utility/modules/builder.hpp>
#include <halmd/utility/modules/concept.hpp>

namespace halmd
{
namespace modules
{

// forward declaration
template <typename T, typename Factory, typename Enable = void>
struct typed_parser_base;

/**
 * This class is constructed upon use of the directive
 * template class module<>, to register concrete (non-abstract)
 * classes with the module mechanism.
 *
 * The module parser iterates through the module's class
 * hierarchy to determine its dependencies and derived-to-base
 * relations. Further it binds any methods needed later during
 * the options parsing stage (options, select) and binds the
 * constructor using a builder.
 */
template <typename T, typename Factory>
struct typed_parser
  : typed_parser_base<T, Factory>
{
    typedef typed_parser_base<T, Factory> Base;
    typedef typename Factory::Registry Registry;
    typedef typename Registry::Graph Graph;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

    /**
     * wrap constructor of the most outer class and attach an
     * untyped base builder shared pointer to the module's vertex
     */
    typed_parser()
    {
        Graph& g = Registry::graph();
        Vertex v = Registry::template vertex<T>();
        property::builder builder(new typed_builder<T, Factory>);
        put(tag::builder(), g, v, builder);
    }
};

/**
 * This class parses a derived class in a module's class hierarchy.
 *
 * The template enables itself if T has a Base typedef.
 */
template <typename T, typename Factory>
struct typed_parser_base<T, Factory, typename boost::enable_if<boost::is_object<typename T::_Base> >::type>
  : typed_parser_base<typename T::_Base, Factory>
{
    typedef typename T::_Base BaseT;
    typedef typed_parser_base<BaseT, Factory> Base;
    typedef typename Factory::Registry Registry;

    /**
     * register derived-to-base relation and module dependencies
     * as edges in the dependency graph and bind methods for use
     * in options parsing stage
     */
    typed_parser_base()
    {
        boost::function_requires<ModuleConcept<T> >();
        Registry::template edge<BaseT, T>(property::is_base_of);
        T::depends();
    }
};

/**
 * This class parses the base in a module's class hierarchy.
 *
 * The template enables itself if T does not have a Base typedef.
 */
template <typename T, typename Factory, typename Enable>
struct typed_parser_base
{
    typedef typename Factory::Registry Registry;

    /**
     * register module dependencies as edges in the dependency
     * graph and bind methods for use in options parsing stage
     */
    typed_parser_base()
    {
        boost::function_requires<ModuleConcept<T> >();
        T::depends();
    }
};

} // namespace modules

} // namespace halmd

#endif /* ! HALMD_UTILITY_MODULES_PARSER_HPP */
