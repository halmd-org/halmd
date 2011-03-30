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

#define BOOST_TEST_MODULE dependency_injection
#include <boost/test/unit_test.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/modules.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <test/tools/lua.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/integrators/verlet.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/tools/cuda.hpp>
#endif

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
namespace test_dummy
{
    struct particle
    {
        virtual ~particle() {}
    };

    struct particle_host : particle {};

    struct integrator
    {
        virtual ~integrator() {}
    };

    struct verlet_host
      : integrator
    {
        verlet_host(shared_ptr<particle_host> particle)
          : particle(particle) {}

        shared_ptr<particle_host> particle;
    };

} // namespace test_dummy

/**
 * Test dependency injection using test dummy modules.
 */
BOOST_FIXTURE_TEST_CASE( dummy, lua_test_fixture )
{
    using namespace test_dummy;
    using namespace luabind;

    LUA_REQUIRE( "assert2 = function(...) assert(...) return assert(...) end" );

    module(L)
    [
        class_<particle, shared_ptr<particle> >("particle")
      , class_<particle_host, shared_ptr<particle>, particle>("particle_host")
            .def(constructor<>())
      , class_<integrator, shared_ptr<integrator> >("integrator")
      , class_<verlet_host, shared_ptr<integrator>, integrator>("verlet_host")
            .def(constructor<shared_ptr<particle_host> >())
            .def_readonly("particle", &verlet_host::particle)
    ];

    LUA_CHECK( "particle = assert2(particle_host())" );
    LUA_CHECK( "integrator = assert2(verlet_host(particle))" );
    LUA_CHECK( "assert(integrator.particle)" );
}

/**
 * Test dependency injection using manually registered HALMD modules.
 *
 * This version uses the host modules.
 */
BOOST_FIXTURE_TEST_CASE( libhalmd_mdsim_host, lua_test_fixture )
{
    luaopen_libhalmd_mdsim_particle(L);
    luaopen_libhalmd_mdsim_host_particle(L);
    luaopen_libhalmd_mdsim_box(L);
    luaopen_libhalmd_mdsim_integrator(L);
    luaopen_libhalmd_mdsim_host_integrators_verlet(L);

    LUA_REQUIRE( "assert2 = function(...) assert(...) return assert(...) end" );

    LUA_CHECK( "particle = assert2(libhalmd.mdsim.host.particle_3_({ 1000 }))" );
    LUA_CHECK( "box = assert2(libhalmd.mdsim.box_3_(particle, { 10, 10, 10 }))" );
    LUA_CHECK( "integrator = assert2(libhalmd.mdsim.host.integrators.verlet_3_(particle, box, 0.001))" );
}

#ifdef WITH_CUDA

struct lua_cuda_fixture : lua_test_fixture, set_cuda_device {};

/**
 * Test dependency injection using manually registered HALMD modules.
 *
 * This version uses the GPU modules.
 */
BOOST_FIXTURE_TEST_CASE( libhalmd_mdsim_gpu, lua_cuda_fixture )
{
    luaopen_libhalmd_mdsim_particle(L);
    luaopen_libhalmd_mdsim_gpu_particle(L);
    luaopen_libhalmd_mdsim_box(L);
    luaopen_libhalmd_mdsim_integrator(L);
    luaopen_libhalmd_mdsim_gpu_integrators_verlet(L);
    luaopen_libhalmd_utility_gpu_device(L);

    LUA_REQUIRE( "assert2 = function(...) assert(...) return assert(...) end" );

    LUA_CHECK( "device = assert2(libhalmd.utility.gpu.device({}, 128))" );
    LUA_CHECK( "particle = assert2(libhalmd.mdsim.gpu.particle_3_(device, { 1000 }))" );
    LUA_CHECK( "box = assert2(libhalmd.mdsim.box_3_(particle, { 10, 10, 10 }))" );
    LUA_CHECK( "integrator = assert2(libhalmd.mdsim.gpu.integrators.verlet_3_(particle, box, 0.001))" );
}

#endif /* WITH_CUDA */
