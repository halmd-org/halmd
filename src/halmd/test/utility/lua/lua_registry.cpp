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

#define BOOST_TEST_MODULE lua_registry
#include <boost/test/unit_test.hpp>

#include <boost/mpl/int.hpp>
#include <halmd/utility/lua/lua_registry.hpp>

using namespace halmd;

// work around incompatibility of macros with multiple template parameters
typedef lua_registry<host_, dimension_<2> > lua_registry_host_dimension_2;
typedef lua_registry<dimension_<2>, host_> lua_registry_dimension_2_host;
typedef lua_registry<gpu_, dimension_<2> > lua_registry_gpu_dimension_2;
typedef lua_registry<dimension_<2>, gpu_> lua_registry_dimension_2_gpu;
typedef lua_registry<host_, dimension_<3> > lua_registry_host_dimension_3;
typedef lua_registry<dimension_<3>, host_> lua_registry_dimension_3_host;
typedef lua_registry<gpu_, dimension_<3> > lua_registry_gpu_dimension_3;
typedef lua_registry<dimension_<3>, gpu_> lua_registry_dimension_3_gpu;

BOOST_AUTO_TEST_CASE( compare_registry_memory_addresses )
{
    BOOST_CHECK_EQUAL( lua_registry<>::get(), lua_registry<>::get() );

    BOOST_CHECK_EQUAL( lua_registry<dimension_<2> >::get(), lua_registry<dimension_<2> >::get() );
    BOOST_CHECK_EQUAL( lua_registry<dimension_<3> >::get(), lua_registry<dimension_<3> >::get() );

    BOOST_CHECK_EQUAL( lua_registry<backend_<host_> >::get(), lua_registry<backend_<host_> >::get() );
    BOOST_CHECK_EQUAL( lua_registry<backend_<gpu_> >::get(), lua_registry<backend_<gpu_> >::get() );
    BOOST_CHECK_EQUAL( lua_registry<gpu_>::get(), lua_registry<backend_<gpu_> >::get() );
    BOOST_CHECK_EQUAL( lua_registry<host_>::get(), lua_registry<backend_<host_> >::get() );

    BOOST_CHECK_NE( lua_registry<>::get(), lua_registry<host_>::get() );
    BOOST_CHECK_NE( lua_registry<dimension_<2> >::get(), lua_registry<dimension_<3> >::get() );
    BOOST_CHECK_NE( lua_registry<gpu_>::get(), lua_registry<host_>::get() );
    BOOST_CHECK_NE( lua_registry<dimension_<2> >::get(), lua_registry<boost::mpl::int_<2> >::get() );
    BOOST_CHECK_NE( lua_registry<dimension_<3> >::get(), lua_registry<boost::mpl::int_<3> >::get() );

    BOOST_CHECK_NE( lua_registry<host_>::get(), lua_registry<dimension_<2> >::get() );
    BOOST_CHECK_NE( lua_registry<host_>::get(), lua_registry<dimension_<3> >::get() );
    BOOST_CHECK_NE( lua_registry<gpu_>::get(), lua_registry<dimension_<2> >::get() );
    BOOST_CHECK_NE( lua_registry<gpu_>::get(), lua_registry<dimension_<3> >::get() );

    BOOST_CHECK_EQUAL( lua_registry_host_dimension_2::get(), lua_registry_host_dimension_2::get() );
    BOOST_CHECK_EQUAL( lua_registry_host_dimension_2::get(), lua_registry_dimension_2_host::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_2::get(), lua_registry_gpu_dimension_2::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_2::get(), lua_registry_dimension_2_gpu::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_2::get(), lua_registry_host_dimension_3::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_2::get(), lua_registry_dimension_3_host::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_2::get(), lua_registry_gpu_dimension_3::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_2::get(), lua_registry_dimension_3_gpu::get() );

    BOOST_CHECK_NE( lua_registry_gpu_dimension_2::get(), lua_registry_host_dimension_2::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_2::get(), lua_registry_dimension_2_host::get() );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_2::get(), lua_registry_gpu_dimension_2::get() );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_2::get(), lua_registry_dimension_2_gpu::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_2::get(), lua_registry_host_dimension_3::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_2::get(), lua_registry_dimension_3_host::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_2::get(), lua_registry_gpu_dimension_3::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_2::get(), lua_registry_dimension_3_gpu::get() );

    BOOST_CHECK_NE( lua_registry_host_dimension_3::get(), lua_registry_host_dimension_2::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_3::get(), lua_registry_dimension_2_host::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_3::get(), lua_registry_gpu_dimension_2::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_3::get(), lua_registry_dimension_2_gpu::get() );
    BOOST_CHECK_EQUAL( lua_registry_host_dimension_3::get(), lua_registry_host_dimension_3::get() );
    BOOST_CHECK_EQUAL( lua_registry_host_dimension_3::get(), lua_registry_dimension_3_host::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_3::get(), lua_registry_gpu_dimension_3::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_3::get(), lua_registry_dimension_3_gpu::get() );

    BOOST_CHECK_NE( lua_registry_gpu_dimension_3::get(), lua_registry_host_dimension_2::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_3::get(), lua_registry_dimension_2_host::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_3::get(), lua_registry_gpu_dimension_2::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_3::get(), lua_registry_dimension_2_gpu::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_3::get(), lua_registry_host_dimension_3::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_3::get(), lua_registry_dimension_3_host::get() );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_3::get(), lua_registry_gpu_dimension_3::get() );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_3::get(), lua_registry_dimension_3_gpu::get() );

    BOOST_CHECK_NE( lua_registry_host_dimension_2::get(), lua_registry<host_>::get() );
    BOOST_CHECK_NE( lua_registry_host_dimension_3::get(), lua_registry<host_>::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_2::get(), lua_registry<gpu_>::get() );
    BOOST_CHECK_NE( lua_registry_gpu_dimension_3::get(), lua_registry<gpu_>::get() );

    BOOST_CHECK_NE( lua_registry_dimension_2_host::get(), lua_registry<dimension_<2> >::get() );
    BOOST_CHECK_NE( lua_registry_dimension_3_host::get(), lua_registry<dimension_<3> >::get() );
    BOOST_CHECK_NE( lua_registry_dimension_2_gpu::get(), lua_registry<dimension_<2> >::get() );
    BOOST_CHECK_NE( lua_registry_dimension_3_gpu::get(), lua_registry<dimension_<3> >::get() );
}

static void register_lua(lua_State* L) {}

static void register_lua_host(lua_State* L) {}

static void register_lua_gpu(lua_State* L) {}

template <int dimension>
static void register_lua_dim(lua_State* L) {}

template <int dimension>
static void register_lua_dim_host(lua_State* L) {}

template <int dimension>
static void register_lua_dim_gpu(lua_State* L) {}

static lua_registry<>::iterator dummy1 = (
    lua_registry<>::get()->push_back( &register_lua )
  , lua_registry<host_ >::get()->push_back( &register_lua_host )
  , lua_registry<gpu_ >::get()->push_back( &register_lua_gpu )
  , lua_registry<dimension_<2> >::get()->push_back( &register_lua_dim<2> )
  , lua_registry<dimension_<3> >::get()->push_back( &register_lua_dim<3> )
  , lua_registry<>::get()->end()
);

BOOST_AUTO_TEST_CASE( register_functions )
{
    BOOST_CHECK_EQUAL( lua_registry<>::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<host_>::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<gpu_>::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<dimension_<2> >::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<dimension_<3> >::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry_host_dimension_2::get()->size(), 0u );
    BOOST_CHECK_EQUAL( lua_registry_host_dimension_2::get()->size(), 0u );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_3::get()->size(), 0u );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_3::get()->size(), 0u );

    BOOST_CHECK_EQUAL( lua_registry<>::get()->front(), &register_lua );
    BOOST_CHECK_EQUAL( lua_registry<host_>::get()->front(), &register_lua_host );
    BOOST_CHECK_EQUAL( lua_registry<gpu_>::get()->front(), &register_lua_gpu );
    BOOST_CHECK_EQUAL( lua_registry<>::get()->front(), &register_lua );
    BOOST_CHECK_EQUAL( lua_registry<host_>::get()->front(), &register_lua_host );
    BOOST_CHECK_EQUAL( lua_registry<gpu_>::get()->front(), &register_lua_gpu );

    lua_registry<dimension_<2>, host_ >::get()->push_back( &register_lua_dim_host<2> );
    lua_registry<dimension_<2>, gpu_ >::get()->push_back( &register_lua_dim_gpu<2> );
    lua_registry<dimension_<3>, host_ >::get()->push_back( &register_lua_dim_host<3> );
    lua_registry<dimension_<3>, gpu_ >::get()->push_back( &register_lua_dim_gpu<3> );

    BOOST_CHECK_EQUAL( lua_registry<>::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<host_>::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<gpu_>::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<dimension_<2> >::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry<dimension_<3> >::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry_host_dimension_2::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry_host_dimension_2::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_3::get()->size(), 1u );
    BOOST_CHECK_EQUAL( lua_registry_gpu_dimension_3::get()->size(), 1u );

    for (size_t i = 0; i < 41; ++i) {
        lua_registry<host_, dimension_<2> >::get()->push_back( &register_lua_dim_host<2> );
        lua_registry<gpu_, dimension_<2> >::get()->push_back( &register_lua_dim_gpu<2> );
        lua_registry<host_, dimension_<3> >::get()->push_back( &register_lua_dim_host<3> );
        lua_registry<gpu_, dimension_<3> >::get()->push_back( &register_lua_dim_gpu<3> );
    }

    BOOST_CHECK_EQUAL( std::count(
        lua_registry_host_dimension_2::get()->begin()
      , lua_registry_host_dimension_2::get()->end()
      , &register_lua_dim_host<2>), 42u );
    BOOST_CHECK_EQUAL( std::count(
        lua_registry_gpu_dimension_2::get()->begin()
      , lua_registry_gpu_dimension_2::get()->end()
      , &register_lua_dim_gpu<2>), 42u );
    BOOST_CHECK_EQUAL( std::count(
        lua_registry_host_dimension_3::get()->begin()
      , lua_registry_host_dimension_3::get()->end()
      , &register_lua_dim_host<3>), 42u );
    BOOST_CHECK_EQUAL( std::count(
        lua_registry_gpu_dimension_3::get()->begin()
      , lua_registry_gpu_dimension_3::get()->end()
      , &register_lua_dim_gpu<3>), 42u );
}
