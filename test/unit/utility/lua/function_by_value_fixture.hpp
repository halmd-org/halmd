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

#ifndef TEST_UNIT_UTILITY_LUA_BY_VALUE_FIXTURE_HPP
#define TEST_UNIT_UTILITY_LUA_BY_VALUE_FIXTURE_HPP

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/test/output_test_stream.hpp>

#include <halmd/utility/lua/function.hpp>
#include <test/tools/lua.hpp>

struct function_by_value_fixture : lua_test_fixture
{
    static boost::test_tools::output_test_stream output;

    function_by_value_fixture()
    {
        using namespace luabind;
        module(L, "test")
        [
            def("slot_by_value0", &wrap_slot_by_value0)
          , def("slot_by_value1", &wrap_slot_by_value1)
          , def("slot_by_value2", &wrap_slot_by_value2)
          , def("slot_by_value3", &wrap_slot_by_value3)
          , def("slot_by_value4", &wrap_slot_by_value4)
          , def("slot_by_value5", &wrap_slot_by_value5)
          , def("slot_by_value6", &wrap_slot_by_value6)
          , def("slot_by_value7", &wrap_slot_by_value7)
          , def("slot_by_value8", &wrap_slot_by_value8)
          , def("slot_by_value9", &wrap_slot_by_value9)
          , def("call", &call_slot_by_value0)
          , def("call", &call_slot_by_value1)
          , def("call", &call_slot_by_value2)
          , def("call", &call_slot_by_value3)
          , def("call", &call_slot_by_value4)
          , def("call", &call_slot_by_value5)
          , def("call", &call_slot_by_value6)
          , def("call", &call_slot_by_value7)
          , def("call", &call_slot_by_value8)
          , def("call", &call_slot_by_value9)
        ];
    }

    static int slot_by_value0()
    {
        output << boost::make_tuple(std::string("slot_by_value0"));
        return 1000;
    };

    static int slot_by_value1(int arg1)
    {
        output << boost::make_tuple(std::string("slot_by_value1"), arg1);
        return 1001;
    };

    static int slot_by_value2(int arg1, int arg2)
    {
        output << boost::make_tuple(std::string("slot_by_value2"), arg1, arg2);
        return 1002;
    };

    static int slot_by_value3(int arg1, int arg2, int arg3)
    {
        output << boost::make_tuple(std::string("slot_by_value3"), arg1, arg2, arg3);
        return 1003;
    };

    static int slot_by_value4(int arg1, int arg2, int arg3, int arg4)
    {
        output << boost::make_tuple(std::string("slot_by_value4"), arg1, arg2, arg3, arg4);
        return 1004;
    };

    static int slot_by_value5(int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        output << boost::make_tuple(std::string("slot_by_value5"), arg1, arg2, arg3, arg4, arg5);
        return 1005;
    };

    static int slot_by_value6(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        output << boost::make_tuple(std::string("slot_by_value6"), arg1, arg2, arg3, arg4, arg5, arg6);
        return 1006;
    };

    static int slot_by_value7(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        output << boost::make_tuple(std::string("slot_by_value7"), arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        return 1007;
    };

    static int slot_by_value8(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        output << boost::make_tuple(std::string("slot_by_value8"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        return 1008;
    };

    static int slot_by_value9(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        output << boost::make_tuple(std::string("slot_by_value9"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        return 1009;
    };

    static std::function<int ()> wrap_slot_by_value0()
    {
        return boost::bind(&slot_by_value0);
    }

    static std::function<int (int)> wrap_slot_by_value1()
    {
        return boost::bind(&slot_by_value1, _1);
    }

    static std::function<int (int, int)> wrap_slot_by_value2()
    {
        return boost::bind(&slot_by_value2, _1, _2);
    }

    static std::function<int (int, int, int)> wrap_slot_by_value3()
    {
        return boost::bind(&slot_by_value3, _1, _2, _3);
    }

    static std::function<int (int, int, int, int)> wrap_slot_by_value4()
    {
        return boost::bind(&slot_by_value4, _1, _2, _3, _4);
    }

    static std::function<int (int, int, int, int, int)> wrap_slot_by_value5()
    {
        return boost::bind(&slot_by_value5, _1, _2, _3, _4, _5);
    }

    static std::function<int (int, int, int, int, int, int)> wrap_slot_by_value6()
    {
        return boost::bind(&slot_by_value6, _1, _2, _3, _4, _5, _6);
    }

    static std::function<int (int, int, int, int, int, int, int)> wrap_slot_by_value7()
    {
        return boost::bind(&slot_by_value7, _1, _2, _3, _4, _5, _6, _7);
    }

    static std::function<int (int, int, int, int, int, int, int, int)> wrap_slot_by_value8()
    {
        return boost::bind(&slot_by_value8, _1, _2, _3, _4, _5, _6, _7, _8);
    }

    static std::function<int (int, int, int, int, int, int, int, int, int)> wrap_slot_by_value9()
    {
        return boost::bind(&slot_by_value9, _1, _2, _3, _4, _5, _6, _7, _8, _9);
    }

    static int call_slot_by_value0(std::function<int ()> const& slot_by_value)
    {
        int value = slot_by_value();
        return value;
    }

    static int call_slot_by_value1(std::function<int (int)> const& slot_by_value, int arg1)
    {
        int value = slot_by_value(arg1);
        return value + arg1;
    }

    static int call_slot_by_value2(std::function<int (int, int)> const& slot_by_value, int arg1, int arg2)
    {
        int value = slot_by_value(arg1, arg2);
        return value + arg1 + arg2;
    }

    static int call_slot_by_value3(std::function<int (int, int, int)> const& slot_by_value, int arg1, int arg2, int arg3)
    {
        int value = slot_by_value(arg1, arg2, arg3);
        return value + arg1 + arg2 + arg3;
    }

    static int call_slot_by_value4(std::function<int (int, int, int, int)> const& slot_by_value, int arg1, int arg2, int arg3, int arg4)
    {
        int value = slot_by_value(arg1, arg2, arg3, arg4);
        return value + arg1 + arg2 + arg3 + arg4;
    }

    static int call_slot_by_value5(std::function<int (int, int, int, int, int)> const& slot_by_value, int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        int value = slot_by_value(arg1, arg2, arg3, arg4, arg5);
        return value + arg1 + arg2 + arg3 + arg4 + arg5;
    }

    static int call_slot_by_value6(std::function<int (int, int, int, int, int, int)> const& slot_by_value, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        int value = slot_by_value(arg1, arg2, arg3, arg4, arg5, arg6);
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6;
    }

    static int call_slot_by_value7(std::function<int (int, int, int, int, int, int, int)> const& slot_by_value, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        int value = slot_by_value(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7;
    }

    static int call_slot_by_value8(std::function<int (int, int, int, int, int, int, int, int)> const& slot_by_value, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        int value = slot_by_value(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8;
    }

    static int call_slot_by_value9(std::function<int (int, int, int, int, int, int, int, int, int)> const& slot_by_value, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        int value = slot_by_value(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8 + arg9;
    }
};

boost::test_tools::output_test_stream function_by_value_fixture::output;

#endif /* ! TEST_UNIT_UTILITY_LUA_BY_VALUE_FIXTURE_HPP */
