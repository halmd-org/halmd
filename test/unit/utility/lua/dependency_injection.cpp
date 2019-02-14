/*
 * Copyright © 2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE dependency_injection
#include <boost/test/unit_test.hpp>

#include <halmd/utility/lua/lua.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/lua.hpp>

/**
 * This test checks the dependency injection of C++ modules from Lua.
 *
 * Modules are successively instantiated in Lua using shared pointers,
 * where dependencies of a module are passed to its constructor as
 * arguments. For dependencies which are derived classes, we check that
 * Luabind down- and up-casts to the type of the constructor argument.
 */

using namespace boost;
using namespace std;

/**
 * A hierarchy of test dummy modules.
 */
namespace test_dummy {

struct particle
{
    virtual ~particle() {}
};

struct particle_host : particle {};

struct box
{
    box(std::shared_ptr<test_dummy::particle> particle)
      : particle(particle) {}

    std::shared_ptr<test_dummy::particle> particle;
};

struct integrator
{
    virtual ~integrator() {}
};

struct verlet_host
  : integrator
{
    verlet_host(std::shared_ptr<test_dummy::particle_host> particle)
      : particle(particle) {}

    std::shared_ptr<test_dummy::particle_host> particle;
};

} // namespace test_dummy

struct bind_test_dummy : lua_test_fixture
{
    bind_test_dummy()
    {
        using namespace luaponte;
        using namespace test_dummy;
        module(L, "test_dummy")
        [
            class_<particle>("particle")

          , class_<particle_host, std::shared_ptr<particle>, particle>("particle_host")
                .def(constructor<>())

          , class_<box, std::shared_ptr<box> >("box")
                .def(constructor<std::shared_ptr<particle> >())
                .def_readonly("particle", &box::particle)

          , class_<integrator>("integrator")

          , class_<verlet_host, std::shared_ptr<integrator>, integrator>("verlet_host")
                .def(constructor<std::shared_ptr<particle_host> >())
                .def_readonly("particle", &verlet_host::particle)
        ];
    }
};

/**
 * Test dependency injection using test dummy modules.
 */
BOOST_FIXTURE_TEST_CASE( shared_ptr_casts, bind_test_dummy )
{
    LUA_REQUIRE( "assert2 = function(...) assert(...) return assert(...) end" );

    // create instance of derived class
    LUA_CHECK( "particle = assert2(test_dummy.particle_host())" );
    // pass to constructor as base class shared_ptr
    LUA_CHECK( "box = assert2(test_dummy.box(particle))" );
    LUA_CHECK( "assert(box.particle)" );
    // pass to constructor as dervived class shared_ptr
    LUA_CHECK( "integrator = assert2(test_dummy.verlet_host(particle))" );
    LUA_CHECK( "assert(integrator.particle)" );
}
