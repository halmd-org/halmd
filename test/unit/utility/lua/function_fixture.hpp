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

#ifndef TEST_UNIT_UTILITY_LUA_FIXTURE_HPP
#define TEST_UNIT_UTILITY_LUA_FIXTURE_HPP

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/test/output_test_stream.hpp>

#include <halmd/utility/lua/function.hpp>
#include <test/tools/lua.hpp>

struct function_fixture : lua_test_fixture
{
    static boost::test_tools::output_test_stream output;

    function_fixture()
    {
        using namespace luaponte;
        module(L, "test")
        [
            def("slot0", &wrap_slot0)
          , def("slot1", &wrap_slot1)
          , def("slot2", &wrap_slot2)
          , def("slot3", &wrap_slot3)
          , def("slot4", &wrap_slot4)
          , def("slot5", &wrap_slot5)
          , def("slot6", &wrap_slot6)
          , def("slot7", &wrap_slot7)
          , def("slot8", &wrap_slot8)
          , def("slot9", &wrap_slot9)
          , def("call", &call_slot0)
          , def("call", &call_slot1)
          , def("call", &call_slot2)
          , def("call", &call_slot3)
          , def("call", &call_slot4)
          , def("call", &call_slot5)
          , def("call", &call_slot6)
          , def("call", &call_slot7)
          , def("call", &call_slot8)
          , def("call", &call_slot9)
        ];
    }

    static void slot0()
    {
        output << boost::make_tuple(std::string("slot0"));
    };

    static void slot1(int arg1)
    {
        output << boost::make_tuple(std::string("slot1"), arg1);
    };

    static void slot2(int arg1, int arg2)
    {
        output << boost::make_tuple(std::string("slot2"), arg1, arg2);
    };

    static void slot3(int arg1, int arg2, int arg3)
    {
        output << boost::make_tuple(std::string("slot3"), arg1, arg2, arg3);
    };

    static void slot4(int arg1, int arg2, int arg3, int arg4)
    {
        output << boost::make_tuple(std::string("slot4"), arg1, arg2, arg3, arg4);
    };

    static void slot5(int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        output << boost::make_tuple(std::string("slot5"), arg1, arg2, arg3, arg4, arg5);
    };

    static void slot6(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        output << boost::make_tuple(std::string("slot6"), arg1, arg2, arg3, arg4, arg5, arg6);
    };

    static void slot7(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        output << boost::make_tuple(std::string("slot7"), arg1, arg2, arg3, arg4, arg5, arg6, arg7);
    };

    static void slot8(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        output << boost::make_tuple(std::string("slot8"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    };

    static void slot9(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        output << boost::make_tuple(std::string("slot9"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    };

    static std::function<void ()> wrap_slot0()
    {
        return boost::bind(&slot0);
    }

    static std::function<void (int)> wrap_slot1()
    {
        return boost::bind(&slot1, _1);
    }

    static std::function<void (int, int)> wrap_slot2()
    {
        return boost::bind(&slot2, _1, _2);
    }

    static std::function<void (int, int, int)> wrap_slot3()
    {
        return boost::bind(&slot3, _1, _2, _3);
    }

    static std::function<void (int, int, int, int)> wrap_slot4()
    {
        return boost::bind(&slot4, _1, _2, _3, _4);
    }

    static std::function<void (int, int, int, int, int)> wrap_slot5()
    {
        return boost::bind(&slot5, _1, _2, _3, _4, _5);
    }

    static std::function<void (int, int, int, int, int, int)> wrap_slot6()
    {
        return boost::bind(&slot6, _1, _2, _3, _4, _5, _6);
    }

    static std::function<void (int, int, int, int, int, int, int)> wrap_slot7()
    {
        return boost::bind(&slot7, _1, _2, _3, _4, _5, _6, _7);
    }

    static std::function<void (int, int, int, int, int, int, int, int)> wrap_slot8()
    {
        return boost::bind(&slot8, _1, _2, _3, _4, _5, _6, _7, _8);
    }

    static std::function<void (int, int, int, int, int, int, int, int, int)> wrap_slot9()
    {
        return boost::bind(&slot9, _1, _2, _3, _4, _5, _6, _7, _8, _9);
    }

    static int call_slot0(std::function<void ()> const& slot)
    {
        slot();
        return 0;
    }

    static int call_slot1(std::function<void (int)> const& slot, int arg1)
    {
        slot(arg1);
        return arg1;
    }

    static int call_slot2(std::function<void (int, int)> const& slot, int arg1, int arg2)
    {
        slot(arg1, arg2);
        return arg1 + arg2;
    }

    static int call_slot3(std::function<void (int, int, int)> const& slot, int arg1, int arg2, int arg3)
    {
        slot(arg1, arg2, arg3);
        return arg1 + arg2 + arg3;
    }

    static int call_slot4(std::function<void (int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4)
    {
        slot(arg1, arg2, arg3, arg4);
        return arg1 + arg2 + arg3 + arg4;
    }

    static int call_slot5(std::function<void (int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        slot(arg1, arg2, arg3, arg4, arg5);
        return arg1 + arg2 + arg3 + arg4 + arg5;
    }

    static int call_slot6(std::function<void (int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        slot(arg1, arg2, arg3, arg4, arg5, arg6);
        return arg1 + arg2 + arg3 + arg4 + arg5 + arg6;
    }

    static int call_slot7(std::function<void (int, int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        slot(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        return arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7;
    }

    static int call_slot8(std::function<void (int, int, int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        slot(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        return arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8;
    }

    static int call_slot9(std::function<void (int, int, int, int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        slot(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        return arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8 + arg9;
    }
};

boost::test_tools::output_test_stream function_fixture::output;

#endif /* TEST_UNIT_UTILITY_LUA_FIXTURE_HPP */
