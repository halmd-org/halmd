/*
 * Copyright © 2012  Peter Colberg
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

#define BOOST_TEST_MODULE function
#include <boost/test/unit_test.hpp>

#include <boost/lexical_cast.hpp>

#include <test/tools/ctest.hpp>
#include <test/unit/utility/lua/function_exception_fixture.hpp>
#include <test/unit/utility/lua/function_by_const_ref_fixture.hpp>
#include <test/unit/utility/lua/function_by_ref_fixture.hpp>
#include <test/unit/utility/lua/function_by_value_fixture.hpp>
#include <test/unit/utility/lua/function_fixture.hpp>

using namespace boost;
using namespace std;

BOOST_FIXTURE_TEST_CASE( cpp_function, function_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "assert(  0 == test.call(test.slot0()))" );
        BOOST_CHECK( output.is_equal("(slot0)") );
        LUA_CHECK( "assert( 10 == test.call(test.slot1(), 10))" );
        BOOST_CHECK( output.is_equal("(slot1 10)") );
        LUA_CHECK( "assert( 30 == test.call(test.slot2(), 10, 20))" );
        BOOST_CHECK( output.is_equal("(slot2 10 20)") );
        LUA_CHECK( "assert( 60 == test.call(test.slot3(), 10, 20, 30))" );
        BOOST_CHECK( output.is_equal("(slot3 10 20 30)") );
        LUA_CHECK( "assert(100 == test.call(test.slot4(), 10, 20, 30, 40))" );
        BOOST_CHECK( output.is_equal("(slot4 10 20 30 40)") );
        LUA_CHECK( "assert(150 == test.call(test.slot5(), 10, 20, 30, 40, 50))" );
        BOOST_CHECK( output.is_equal("(slot5 10 20 30 40 50)") );
        LUA_CHECK( "assert(210 == test.call(test.slot6(), 10, 20, 30, 40, 50, 60))" );
        BOOST_CHECK( output.is_equal("(slot6 10 20 30 40 50 60)") );
        LUA_CHECK( "assert(280 == test.call(test.slot7(), 10, 20, 30, 40, 50, 60, 70))" );
        BOOST_CHECK( output.is_equal("(slot7 10 20 30 40 50 60 70)") );
        LUA_CHECK( "assert(360 == test.call(test.slot8(), 10, 20, 30, 40, 50, 60, 70, 80))" );
        BOOST_CHECK( output.is_equal("(slot8 10 20 30 40 50 60 70 80)") );
        LUA_CHECK( "assert(450 == test.call(test.slot9(), 10, 20, 30, 40, 50, 60, 70, 80, 90))" );
        BOOST_CHECK( output.is_equal("(slot9 10 20 30 40 50 60 70 80 90)") );
    }
}

BOOST_FIXTURE_TEST_CASE( cpp_function_by_value, function_by_value_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "assert(1000 == test.call(test.slot_by_value0()))" );
        BOOST_CHECK( output.is_equal("(slot_by_value0)") );
        LUA_CHECK( "assert(1011 == test.call(test.slot_by_value1(), 10))" );
        BOOST_CHECK( output.is_equal("(slot_by_value1 10)") );
        LUA_CHECK( "assert(1032 == test.call(test.slot_by_value2(), 10, 20))" );
        BOOST_CHECK( output.is_equal("(slot_by_value2 10 20)") );
        LUA_CHECK( "assert(1063 == test.call(test.slot_by_value3(), 10, 20, 30))" );
        BOOST_CHECK( output.is_equal("(slot_by_value3 10 20 30)") );
        LUA_CHECK( "assert(1104 == test.call(test.slot_by_value4(), 10, 20, 30, 40))" );
        BOOST_CHECK( output.is_equal("(slot_by_value4 10 20 30 40)") );
        LUA_CHECK( "assert(1155 == test.call(test.slot_by_value5(), 10, 20, 30, 40, 50))" );
        BOOST_CHECK( output.is_equal("(slot_by_value5 10 20 30 40 50)") );
        LUA_CHECK( "assert(1216 == test.call(test.slot_by_value6(), 10, 20, 30, 40, 50, 60))" );
        BOOST_CHECK( output.is_equal("(slot_by_value6 10 20 30 40 50 60)") );
        LUA_CHECK( "assert(1287 == test.call(test.slot_by_value7(), 10, 20, 30, 40, 50, 60, 70))" );
        BOOST_CHECK( output.is_equal("(slot_by_value7 10 20 30 40 50 60 70)") );
        LUA_CHECK( "assert(1368 == test.call(test.slot_by_value8(), 10, 20, 30, 40, 50, 60, 70, 80))" );
        BOOST_CHECK( output.is_equal("(slot_by_value8 10 20 30 40 50 60 70 80)") );
        LUA_CHECK( "assert(1459 == test.call(test.slot_by_value9(), 10, 20, 30, 40, 50, 60, 70, 80, 90))" );
        BOOST_CHECK( output.is_equal("(slot_by_value9 10 20 30 40 50 60 70 80 90)") );
    }
}

BOOST_FIXTURE_TEST_CASE( cpp_function_by_const_ref, function_by_const_ref_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0000 == test.call(test.slot_by_const_ref0()))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref0)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0011 == test.call(test.slot_by_const_ref1(), 10))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref1 10)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0032 == test.call(test.slot_by_const_ref2(), 10, 20))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref2 10 20)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0063 == test.call(test.slot_by_const_ref3(), 10, 20, 30))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref3 10 20 30)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0104 == test.call(test.slot_by_const_ref4(), 10, 20, 30, 40))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref4 10 20 30 40)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0155 == test.call(test.slot_by_const_ref5(), 10, 20, 30, 40, 50))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref5 10 20 30 40 50)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0216 == test.call(test.slot_by_const_ref6(), 10, 20, 30, 40, 50, 60))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref6 10 20 30 40 50 60)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0287 == test.call(test.slot_by_const_ref7(), 10, 20, 30, 40, 50, 60, 70))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref7 10 20 30 40 50 60 70)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0368 == test.call(test.slot_by_const_ref8(), 10, 20, 30, 40, 50, 60, 70, 80))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref8 10 20 30 40 50 60 70 80)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "0459 == test.call(test.slot_by_const_ref9(), 10, 20, 30, 40, 50, 60, 70, 80, 90))" );
        BOOST_CHECK( output.is_equal("(slot_by_const_ref9 10 20 30 40 50 60 70 80 90)") );
    }
}

BOOST_FIXTURE_TEST_CASE( cpp_function_by_ref, function_by_ref_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00000 == test.call(test.slot_by_ref0()))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref0)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00011 == test.call(test.slot_by_ref1(), 10))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref1 10)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00032 == test.call(test.slot_by_ref2(), 10, 20))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref2 10 20)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00063 == test.call(test.slot_by_ref3(), 10, 20, 30))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref3 10 20 30)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00104 == test.call(test.slot_by_ref4(), 10, 20, 30, 40))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref4 10 20 30 40)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00155 == test.call(test.slot_by_ref5(), 10, 20, 30, 40, 50))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref5 10 20 30 40 50)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00216 == test.call(test.slot_by_ref6(), 10, 20, 30, 40, 50, 60))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref6 10 20 30 40 50 60)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00287 == test.call(test.slot_by_ref7(), 10, 20, 30, 40, 50, 60, 70))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref7 10 20 30 40 50 60 70)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00368 == test.call(test.slot_by_ref8(), 10, 20, 30, 40, 50, 60, 70, 80))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref8 10 20 30 40 50 60 70 80)") );
        LUA_CHECK( "assert(" + lexical_cast<string>(i) + "00459 == test.call(test.slot_by_ref9(), 10, 20, 30, 40, 50, 60, 70, 80, 90))" );
        BOOST_CHECK( output.is_equal("(slot_by_ref9 10 20 30 40 50 60 70 80 90)") );
    }
}

BOOST_FIXTURE_TEST_CASE( cpp_function_exception, function_exception_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "test.catch_exception(test.throw_exception0())" );
        BOOST_CHECK( output.is_equal("(exception0)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception1(), 10)" );
        BOOST_CHECK( output.is_equal("(exception1 10)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception2(), 10, 20)" );
        BOOST_CHECK( output.is_equal("(exception2 10 20)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception3(), 10, 20, 30)" );
        BOOST_CHECK( output.is_equal("(exception3 10 20 30)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception4(), 10, 20, 30, 40)" );
        BOOST_CHECK( output.is_equal("(exception4 10 20 30 40)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception5(), 10, 20, 30, 40, 50)" );
        BOOST_CHECK( output.is_equal("(exception5 10 20 30 40 50)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception6(), 10, 20, 30, 40, 50, 60)" );
        BOOST_CHECK( output.is_equal("(exception6 10 20 30 40 50 60)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception7(), 10, 20, 30, 40, 50, 60, 70)" );
        BOOST_CHECK( output.is_equal("(exception7 10 20 30 40 50 60 70)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception8(), 10, 20, 30, 40, 50, 60, 70, 80)" );
        BOOST_CHECK( output.is_equal("(exception8 10 20 30 40 50 60 70 80)") );
        LUA_CHECK( "test.catch_exception(test.throw_exception9(), 10, 20, 30, 40, 50, 60, 70, 80, 90)" );
        BOOST_CHECK( output.is_equal("(exception9 10 20 30 40 50 60 70 80 90)") );
    }
}

BOOST_FIXTURE_TEST_CASE( lua_function_slot, function_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "assert(  0 == test.call(function(...) output = table.concat({'slot0', ...}, ' ') end))" );
        LUA_CHECK( "assert(output == 'slot0')" );
        LUA_CHECK( "assert( 10 == test.call(function(...) output = table.concat({'slot1', ...}, ' ') end, 10))" );
        LUA_CHECK( "assert(output == 'slot1 10')" );
        LUA_CHECK( "assert( 30 == test.call(function(...) output = table.concat({'slot2', ...}, ' ') end, 10, 20))" );
        LUA_CHECK( "assert(output == 'slot2 10 20')" );
        LUA_CHECK( "assert( 60 == test.call(function(...) output = table.concat({'slot3', ...}, ' ') end, 10, 20, 30))" );
        LUA_CHECK( "assert(output == 'slot3 10 20 30')" );
        LUA_CHECK( "assert(100 == test.call(function(...) output = table.concat({'slot4', ...}, ' ') end, 10, 20, 30, 40))" );
        LUA_CHECK( "assert(output == 'slot4 10 20 30 40')" );
        LUA_CHECK( "assert(150 == test.call(function(...) output = table.concat({'slot5', ...}, ' ') end, 10, 20, 30, 40, 50))" );
        LUA_CHECK( "assert(output == 'slot5 10 20 30 40 50')" );
        LUA_CHECK( "assert(210 == test.call(function(...) output = table.concat({'slot6', ...}, ' ') end, 10, 20, 30, 40, 50, 60))" );
        LUA_CHECK( "assert(output == 'slot6 10 20 30 40 50 60')" );
        LUA_CHECK( "assert(280 == test.call(function(...) output = table.concat({'slot7', ...}, ' ') end, 10, 20, 30, 40, 50, 60, 70))" );
        LUA_CHECK( "assert(output == 'slot7 10 20 30 40 50 60 70')" );
        LUA_CHECK( "assert(360 == test.call(function(...) output = table.concat({'slot8', ...}, ' ') end, 10, 20, 30, 40, 50, 60, 70, 80))" );
        LUA_CHECK( "assert(output == 'slot8 10 20 30 40 50 60 70 80')" );
        LUA_CHECK( "assert(450 == test.call(function(...) output = table.concat({'slot9', ...}, ' ') end, 10, 20, 30, 40, 50, 60, 70, 80, 90))" );
        LUA_CHECK( "assert(output == 'slot9 10 20 30 40 50 60 70 80 90')" );
    }
}

BOOST_FIXTURE_TEST_CASE( lua_function_by_value, function_by_value_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "assert(2000 == test.call(function(...) output = table.concat({'slot_by_value0', ...}, ' ') return 2000 end))" );
        LUA_CHECK( "assert(output == 'slot_by_value0')" );
        LUA_CHECK( "assert(2011 == test.call(function(...) output = table.concat({'slot_by_value1', ...}, ' ') return 2001 end, 10))" );
        LUA_CHECK( "assert(output == 'slot_by_value1 10')" );
        LUA_CHECK( "assert(2032 == test.call(function(...) output = table.concat({'slot_by_value2', ...}, ' ') return 2002 end, 10, 20))" );
        LUA_CHECK( "assert(output == 'slot_by_value2 10 20')" );
        LUA_CHECK( "assert(2063 == test.call(function(...) output = table.concat({'slot_by_value3', ...}, ' ') return 2003 end, 10, 20, 30))" );
        LUA_CHECK( "assert(output == 'slot_by_value3 10 20 30')" );
        LUA_CHECK( "assert(2104 == test.call(function(...) output = table.concat({'slot_by_value4', ...}, ' ') return 2004 end, 10, 20, 30, 40))" );
        LUA_CHECK( "assert(output == 'slot_by_value4 10 20 30 40')" );
        LUA_CHECK( "assert(2155 == test.call(function(...) output = table.concat({'slot_by_value5', ...}, ' ') return 2005 end, 10, 20, 30, 40, 50))" );
        LUA_CHECK( "assert(output == 'slot_by_value5 10 20 30 40 50')" );
        LUA_CHECK( "assert(2216 == test.call(function(...) output = table.concat({'slot_by_value6', ...}, ' ') return 2006 end, 10, 20, 30, 40, 50, 60))" );
        LUA_CHECK( "assert(output == 'slot_by_value6 10 20 30 40 50 60')" );
        LUA_CHECK( "assert(2287 == test.call(function(...) output = table.concat({'slot_by_value7', ...}, ' ') return 2007 end, 10, 20, 30, 40, 50, 60, 70))" );
        LUA_CHECK( "assert(output == 'slot_by_value7 10 20 30 40 50 60 70')" );
        LUA_CHECK( "assert(2368 == test.call(function(...) output = table.concat({'slot_by_value8', ...}, ' ') return 2008 end, 10, 20, 30, 40, 50, 60, 70, 80))" );
        LUA_CHECK( "assert(output == 'slot_by_value8 10 20 30 40 50 60 70 80')" );
        LUA_CHECK( "assert(2459 == test.call(function(...) output = table.concat({'slot_by_value9', ...}, ' ') return 2009 end, 10, 20, 30, 40, 50, 60, 70, 80, 90))" );
        LUA_CHECK( "assert(output == 'slot_by_value9 10 20 30 40 50 60 70 80 90')" );
    }
}

BOOST_FIXTURE_TEST_CASE( lua_function_not_by_const_ref, function_by_const_ref_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "status, output = pcall(test.call, function(...) end)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50, 60)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50, 60, 70)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50, 60, 70, 80, 90)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
    }
}

BOOST_FIXTURE_TEST_CASE( lua_function_not_by_ref, function_by_ref_fixture )
{
    for (int i = 1; i <= 3; ++i) {
        BOOST_TEST_MESSAGE( "iteration " << i);
        LUA_CHECK( "status, output = pcall(test.call, function(...) end)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50, 60)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50, 60, 70)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
        LUA_CHECK( "status, output = pcall(test.call, function(...) end, 10, 20, 30, 40, 50, 60, 70, 80, 90)" );
        LUA_CHECK( "assert(output:match('No matching overload found'))" );
    }
}
