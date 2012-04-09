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

#ifndef TEST_UNIT_UTILITY_LUA_BY_CONST_REF_FIXTURE_HPP
#define TEST_UNIT_UTILITY_LUA_BY_CONST_REF_FIXTURE_HPP

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/test/output_test_stream.hpp>

#include <halmd/utility/lua/function.hpp>
#include <test/tools/lua.hpp>

struct function_by_const_ref_fixture : lua_test_fixture
{
    static boost::test_tools::output_test_stream output;

    function_by_const_ref_fixture()
    {
        using namespace luabind;
        module(L, "test")
        [
            def("slot_by_const_ref0", &wrap_slot_by_const_ref0)
          , def("slot_by_const_ref1", &wrap_slot_by_const_ref1)
          , def("slot_by_const_ref2", &wrap_slot_by_const_ref2)
          , def("slot_by_const_ref3", &wrap_slot_by_const_ref3)
          , def("slot_by_const_ref4", &wrap_slot_by_const_ref4)
          , def("slot_by_const_ref5", &wrap_slot_by_const_ref5)
          , def("slot_by_const_ref6", &wrap_slot_by_const_ref6)
          , def("slot_by_const_ref7", &wrap_slot_by_const_ref7)
          , def("slot_by_const_ref8", &wrap_slot_by_const_ref8)
          , def("slot_by_const_ref9", &wrap_slot_by_const_ref9)
          , def("call", &call_slot_by_const_ref0)
          , def("call", &call_slot_by_const_ref1)
          , def("call", &call_slot_by_const_ref2)
          , def("call", &call_slot_by_const_ref3)
          , def("call", &call_slot_by_const_ref4)
          , def("call", &call_slot_by_const_ref5)
          , def("call", &call_slot_by_const_ref6)
          , def("call", &call_slot_by_const_ref7)
          , def("call", &call_slot_by_const_ref8)
          , def("call", &call_slot_by_const_ref9)
        ];
    }

    static int const& slot_by_const_ref0()
    {
        output << boost::make_tuple(std::string("slot_by_const_ref0"));
        static int value = 0;
        return value;
    };

    static int const& slot_by_const_ref1(int arg1)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref1"), arg1);
        static int value = 1;
        return value;
    };

    static int const& slot_by_const_ref2(int arg1, int arg2)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref2"), arg1, arg2);
        static int value = 2;
        return value;
    };

    static int const& slot_by_const_ref3(int arg1, int arg2, int arg3)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref3"), arg1, arg2, arg3);
        static int value = 3;
        return value;
    };

    static int const& slot_by_const_ref4(int arg1, int arg2, int arg3, int arg4)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref4"), arg1, arg2, arg3, arg4);
        static int value = 4;
        return value;
    };

    static int const& slot_by_const_ref5(int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref5"), arg1, arg2, arg3, arg4, arg5);
        static int value = 5;
        return value;
    };

    static int const& slot_by_const_ref6(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref6"), arg1, arg2, arg3, arg4, arg5, arg6);
        static int value = 6;
        return value;
    };

    static int const& slot_by_const_ref7(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref7"), arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        static int value = 7;
        return value;
    };

    static int const& slot_by_const_ref8(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref8"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        static int value = 8;
        return value;
    };

    static int const& slot_by_const_ref9(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        output << boost::make_tuple(std::string("slot_by_const_ref9"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        static int value = 9;
        return value;
    };

    static boost::function<int const& ()> wrap_slot_by_const_ref0()
    {
        return boost::bind(&slot_by_const_ref0);
    }

    static boost::function<int const& (int)> wrap_slot_by_const_ref1()
    {
        return boost::bind(&slot_by_const_ref1, _1);
    }

    static boost::function<int const& (int, int)> wrap_slot_by_const_ref2()
    {
        return boost::bind(&slot_by_const_ref2, _1, _2);
    }

    static boost::function<int const& (int, int, int)> wrap_slot_by_const_ref3()
    {
        return boost::bind(&slot_by_const_ref3, _1, _2, _3);
    }

    static boost::function<int const& (int, int, int, int)> wrap_slot_by_const_ref4()
    {
        return boost::bind(&slot_by_const_ref4, _1, _2, _3, _4);
    }

    static boost::function<int const& (int, int, int, int, int)> wrap_slot_by_const_ref5()
    {
        return boost::bind(&slot_by_const_ref5, _1, _2, _3, _4, _5);
    }

    static boost::function<int const& (int, int, int, int, int, int)> wrap_slot_by_const_ref6()
    {
        return boost::bind(&slot_by_const_ref6, _1, _2, _3, _4, _5, _6);
    }

    static boost::function<int const& (int, int, int, int, int, int, int)> wrap_slot_by_const_ref7()
    {
        return boost::bind(&slot_by_const_ref7, _1, _2, _3, _4, _5, _6, _7);
    }

    static boost::function<int const& (int, int, int, int, int, int, int, int)> wrap_slot_by_const_ref8()
    {
        return boost::bind(&slot_by_const_ref8, _1, _2, _3, _4, _5, _6, _7, _8);
    }

    static boost::function<int const& (int, int, int, int, int, int, int, int, int)> wrap_slot_by_const_ref9()
    {
        return boost::bind(&slot_by_const_ref9, _1, _2, _3, _4, _5, _6, _7, _8, _9);
    }

    static int call_slot_by_const_ref0(boost::function<int const& ()> const& slot_by_const_ref)
    {
        int const& value = slot_by_const_ref();
        const_cast<int&>(value) += 10000;
        return value;
    }

    static int call_slot_by_const_ref1(boost::function<int const& (int)> const& slot_by_const_ref, int arg1)
    {
        int const& value = slot_by_const_ref(arg1);
        const_cast<int&>(value) += 10000;
        return value + arg1;
    }

    static int call_slot_by_const_ref2(boost::function<int const& (int, int)> const& slot_by_const_ref, int arg1, int arg2)
    {
        int const& value = slot_by_const_ref(arg1, arg2);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2;
    }

    static int call_slot_by_const_ref3(boost::function<int const& (int, int, int)> const& slot_by_const_ref, int arg1, int arg2, int arg3)
    {
        int const& value = slot_by_const_ref(arg1, arg2, arg3);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2 + arg3;
    }

    static int call_slot_by_const_ref4(boost::function<int const& (int, int, int, int)> const& slot_by_const_ref, int arg1, int arg2, int arg3, int arg4)
    {
        int const& value = slot_by_const_ref(arg1, arg2, arg3, arg4);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2 + arg3 + arg4;
    }

    static int call_slot_by_const_ref5(boost::function<int const& (int, int, int, int, int)> const& slot_by_const_ref, int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        int const& value = slot_by_const_ref(arg1, arg2, arg3, arg4, arg5);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5;
    }

    static int call_slot_by_const_ref6(boost::function<int const& (int, int, int, int, int, int)> const& slot_by_const_ref, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        int const& value = slot_by_const_ref(arg1, arg2, arg3, arg4, arg5, arg6);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6;
    }

    static int call_slot_by_const_ref7(boost::function<int const& (int, int, int, int, int, int, int)> const& slot_by_const_ref, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        int const& value = slot_by_const_ref(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7;
    }

    static int call_slot_by_const_ref8(boost::function<int const& (int, int, int, int, int, int, int, int)> const& slot_by_const_ref, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        int const& value = slot_by_const_ref(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8;
    }

    static int call_slot_by_const_ref9(boost::function<int const& (int, int, int, int, int, int, int, int, int)> const& slot_by_const_ref, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        int const& value = slot_by_const_ref(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        const_cast<int&>(value) += 10000;
        return value + arg1 + arg2 + arg3 + arg4 + arg5 + arg6 + arg7 + arg8 + arg9;
    }
};

boost::test_tools::output_test_stream function_by_const_ref_fixture::output;

#endif /* ! TEST_UNIT_UTILITY_LUA_BY_CONST_REF_FIXTURE_HPP */
