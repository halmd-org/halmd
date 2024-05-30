/*
 * Copyright © 2012  Peter Colberg
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

#ifndef TEST_UNIT_UTILITY_LUA_BY_REF_FIXTURE_HPP
#define TEST_UNIT_UTILITY_LUA_BY_REF_FIXTURE_HPP

#include <tuple>
#include <boost/bind.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 107100
# include <boost/test/tools/output_test_stream.hpp>
#else
# include <boost/test/output_test_stream.hpp>
#endif

#include <halmd/utility/lua/function.hpp>
#include <test/tools/lua.hpp>
#include <test/tools/tuple_io.hpp>

struct function_by_ref_fixture : lua_test_fixture
{
    static boost::test_tools::output_test_stream output;

    function_by_ref_fixture()
    {
        using namespace luaponte;
        module(L, "test")
        [
            def("slot_by_ref0", &wrap_slot_by_ref0)
          , def("slot_by_ref1", &wrap_slot_by_ref1)
          , def("slot_by_ref2", &wrap_slot_by_ref2)
          , def("slot_by_ref3", &wrap_slot_by_ref3)
          , def("slot_by_ref4", &wrap_slot_by_ref4)
          , def("slot_by_ref5", &wrap_slot_by_ref5)
          , def("slot_by_ref6", &wrap_slot_by_ref6)
          , def("slot_by_ref7", &wrap_slot_by_ref7)
          , def("slot_by_ref8", &wrap_slot_by_ref8)
          , def("call", &call_slot_by_ref0)
          , def("call", &call_slot_by_ref1)
          , def("call", &call_slot_by_ref2)
          , def("call", &call_slot_by_ref3)
          , def("call", &call_slot_by_ref4)
          , def("call", &call_slot_by_ref5)
          , def("call", &call_slot_by_ref6)
          , def("call", &call_slot_by_ref7)
          , def("call", &call_slot_by_ref8)
        ];
    }

    static int& slot_by_ref0()
    {
        output << std::make_tuple(std::string("slot_by_ref0"));
        static int value = 0;
        return value;
    };

    static int& slot_by_ref1(int arg1)
    {
        output << std::make_tuple(std::string("slot_by_ref1"), arg1);
        static int value = 1;
        return value;
    };

    static int& slot_by_ref2(int arg1, int arg2)
    {
        output << std::make_tuple(std::string("slot_by_ref2"), arg1, arg2);
        static int value = 2;
        return value;
    };

    static int& slot_by_ref3(int arg1, int arg2, int arg3)
    {
        output << std::make_tuple(std::string("slot_by_ref3"), arg1, arg2, arg3);
        static int value = 3;
        return value;
    };

    static int& slot_by_ref4(int arg1, int arg2, int arg3, int arg4)
    {
        output << std::make_tuple(std::string("slot_by_ref4"), arg1, arg2, arg3, arg4);
        static int value = 4;
        return value;
    };

    static int& slot_by_ref5(int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        output << std::make_tuple(std::string("slot_by_ref5"), arg1, arg2, arg3, arg4, arg5);
        static int value = 5;
        return value;
    };

    static int& slot_by_ref6(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        output << std::make_tuple(std::string("slot_by_ref6"), arg1, arg2, arg3, arg4, arg5, arg6);
        static int value = 6;
        return value;
    };

    static int& slot_by_ref7(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        output << std::make_tuple(std::string("slot_by_ref7"), arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        static int value = 7;
        return value;
    };

    static int& slot_by_ref8(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        output << std::make_tuple(std::string("slot_by_ref8"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        static int value = 8;
        return value;
    };

    static std::function<int& ()> wrap_slot_by_ref0()
    {
        return boost::bind(&slot_by_ref0);
    }

    static std::function<int& (int)> wrap_slot_by_ref1()
    {
        return boost::bind(&slot_by_ref1, _1);
    }

    static std::function<int& (int, int)> wrap_slot_by_ref2()
    {
        return boost::bind(&slot_by_ref2, _1, _2);
    }

    static std::function<int& (int, int, int)> wrap_slot_by_ref3()
    {
        return boost::bind(&slot_by_ref3, _1, _2, _3);
    }

    static std::function<int& (int, int, int, int)> wrap_slot_by_ref4()
    {
        return boost::bind(&slot_by_ref4, _1, _2, _3, _4);
    }

    static std::function<int& (int, int, int, int, int)> wrap_slot_by_ref5()
    {
        return boost::bind(&slot_by_ref5, _1, _2, _3, _4, _5);
    }

    static std::function<int& (int, int, int, int, int, int)> wrap_slot_by_ref6()
    {
        return boost::bind(&slot_by_ref6, _1, _2, _3, _4, _5, _6);
    }

    static std::function<int& (int, int, int, int, int, int, int)> wrap_slot_by_ref7()
    {
        return boost::bind(&slot_by_ref7, _1, _2, _3, _4, _5, _6, _7);
    }

    static std::function<int& (int, int, int, int, int, int, int, int)> wrap_slot_by_ref8()
    {
        return boost::bind(&slot_by_ref8, _1, _2, _3, _4, _5, _6, _7, _8);
    }

    static int call_slot_by_ref0(std::function<int& ()> const& slot_by_ref)
    {
        int& value = slot_by_ref();
        value += 100000;
        return value;
    }

    static int call_slot_by_ref1(std::function<int& (int)> const& slot_by_ref, int arg1)
    {
        int& value = slot_by_ref(arg1);
        value += 100000;
        return value + arg1;
    }

    static int call_slot_by_ref2(std::function<int& (int, int)> const& slot_by_ref, int arg1, int arg2)
    {
        int& value = slot_by_ref(arg1, arg2);
        value += 100000;
        return value + arg1 + arg2;
    }

    static int call_slot_by_ref3(std::function<int& (int, int, int)> const& slot_by_ref, int arg1, int arg2, int arg3)
    {
        int& value = slot_by_ref(arg1, arg2, arg3);
        value += 100000;
        return value + arg1 + arg2 + arg3;
    }

    static int call_slot_by_ref4(std::function<int& (int, int, int, int)> const& slot_by_ref, int arg1, int arg2, int arg3, int arg4)
    {
        int& value = slot_by_ref(arg1, arg2, arg3, arg4);
        value += 100000;
        return value + arg1 + arg2 + arg3 + arg4;
    }

    static int call_slot_by_ref5(std::function<int& (int, int, int, int, int)> const& slot_by_ref, int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        int& value = slot_by_ref(arg1, arg2, arg3, arg4, arg5);
        value += 100000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5;
    }

    static int call_slot_by_ref6(std::function<int& (int, int, int, int, int, int)> const& slot_by_ref, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        int& value = slot_by_ref(arg1, arg2, arg3, arg4, arg5, arg6);
        value += 100000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6;
    }

    static int call_slot_by_ref7(std::function<int& (int, int, int, int, int, int, int)> const& slot_by_ref, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        int& value = slot_by_ref(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        value += 100000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7;
    }

    static int call_slot_by_ref8(std::function<int& (int, int, int, int, int, int, int, int)> const& slot_by_ref, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        int& value = slot_by_ref(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        value += 100000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8;
    }
};

boost::test_tools::output_test_stream function_by_ref_fixture::output;

#endif /* ! TEST_UNIT_UTILITY_LUA_BY_REF_FIXTURE_HPP */
