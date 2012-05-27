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

#ifndef TEST_UNIT_UTILITY_LUA_EXCEPTION_FIXTURE_HPP
#define TEST_UNIT_UTILITY_LUA_EXCEPTION_FIXTURE_HPP

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/test/output_test_stream.hpp>

#include <halmd/utility/lua/function.hpp>
#include <test/tools/lua.hpp>

struct function_exception_fixture : lua_test_fixture
{
    static boost::test_tools::output_test_stream output;

    function_exception_fixture()
    {
        using namespace luabind;
        module(L, "test")
        [
            def("throw_exception0", &wrap_throw_exception0)
          , def("throw_exception1", &wrap_throw_exception1)
          , def("throw_exception2", &wrap_throw_exception2)
          , def("throw_exception3", &wrap_throw_exception3)
          , def("throw_exception4", &wrap_throw_exception4)
          , def("throw_exception5", &wrap_throw_exception5)
          , def("throw_exception6", &wrap_throw_exception6)
          , def("throw_exception7", &wrap_throw_exception7)
          , def("throw_exception8", &wrap_throw_exception8)
          , def("throw_exception9", &wrap_throw_exception9)
          , def("catch_exception", &catch_exception0)
          , def("catch_exception", &catch_exception1)
          , def("catch_exception", &catch_exception2)
          , def("catch_exception", &catch_exception3)
          , def("catch_exception", &catch_exception4)
          , def("catch_exception", &catch_exception5)
          , def("catch_exception", &catch_exception6)
          , def("catch_exception", &catch_exception7)
          , def("catch_exception", &catch_exception8)
          , def("catch_exception", &catch_exception9)
        ];
    }

    static void throw_exception0()
    {
        throw boost::make_tuple(std::string("exception0"));
    };

    static void throw_exception1(int arg1)
    {
        throw boost::make_tuple(std::string("exception1"), arg1);
    };

    static void throw_exception2(int arg1, int arg2)
    {
        throw boost::make_tuple(std::string("exception2"), arg1, arg2);
    };

    static void throw_exception3(int arg1, int arg2, int arg3)
    {
        throw boost::make_tuple(std::string("exception3"), arg1, arg2, arg3);
    };

    static void throw_exception4(int arg1, int arg2, int arg3, int arg4)
    {
        throw boost::make_tuple(std::string("exception4"), arg1, arg2, arg3, arg4);
    };

    static void throw_exception5(int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        throw boost::make_tuple(std::string("exception5"), arg1, arg2, arg3, arg4, arg5);
    };

    static void throw_exception6(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        throw boost::make_tuple(std::string("exception6"), arg1, arg2, arg3, arg4, arg5, arg6);
    };

    static void throw_exception7(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        throw boost::make_tuple(std::string("exception7"), arg1, arg2, arg3, arg4, arg5, arg6, arg7);
    };

    static void throw_exception8(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        throw boost::make_tuple(std::string("exception8"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
    };

    static void throw_exception9(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        throw boost::make_tuple(std::string("exception9"), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
    };

    static std::function<void ()> wrap_throw_exception0()
    {
        return boost::bind(&throw_exception0);
    }

    static std::function<void (int)> wrap_throw_exception1()
    {
        return boost::bind(&throw_exception1, _1);
    }

    static std::function<void (int, int)> wrap_throw_exception2()
    {
        return boost::bind(&throw_exception2, _1, _2);
    }

    static std::function<void (int, int, int)> wrap_throw_exception3()
    {
        return boost::bind(&throw_exception3, _1, _2, _3);
    }

    static std::function<void (int, int, int, int)> wrap_throw_exception4()
    {
        return boost::bind(&throw_exception4, _1, _2, _3, _4);
    }

    static std::function<void (int, int, int, int, int)> wrap_throw_exception5()
    {
        return boost::bind(&throw_exception5, _1, _2, _3, _4, _5);
    }

    static std::function<void (int, int, int, int, int, int)> wrap_throw_exception6()
    {
        return boost::bind(&throw_exception6, _1, _2, _3, _4, _5, _6);
    }

    static std::function<void (int, int, int, int, int, int, int)> wrap_throw_exception7()
    {
        return boost::bind(&throw_exception7, _1, _2, _3, _4, _5, _6, _7);
    }

    static std::function<void (int, int, int, int, int, int, int, int)> wrap_throw_exception8()
    {
        return boost::bind(&throw_exception8, _1, _2, _3, _4, _5, _6, _7, _8);
    }

    static std::function<void (int, int, int, int, int, int, int, int, int)> wrap_throw_exception9()
    {
        return boost::bind(&throw_exception9, _1, _2, _3, _4, _5, _6, _7, _8, _9);
    }

    static void catch_exception0(std::function<void ()> const& slot)
    {
        try {
            slot();
        }
        catch (boost::tuple<std::string> const& t)
        {
            output << t;
        }
    }

    static void catch_exception1(std::function<void (int)> const& slot, int arg1)
    {
        try {
            slot(arg1);
        }
        catch (boost::tuple<std::string, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception2(std::function<void (int, int)> const& slot, int arg1, int arg2)
    {
        try {
            slot(arg1, arg2);
        }
        catch (boost::tuple<std::string, int, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception3(std::function<void (int, int, int)> const& slot, int arg1, int arg2, int arg3)
    {
        try {
            slot(arg1, arg2, arg3);
        }
        catch (boost::tuple<std::string, int, int, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception4(std::function<void (int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4)
    {
        try {
            slot(arg1, arg2, arg3, arg4);
        }
        catch (boost::tuple<std::string, int, int, int, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception5(std::function<void (int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        try {
            slot(arg1, arg2, arg3, arg4, arg5);
        }
        catch (boost::tuple<std::string, int, int, int, int, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception6(std::function<void (int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        try {
            slot(arg1, arg2, arg3, arg4, arg5, arg6);
        }
        catch (boost::tuple<std::string, int, int, int, int, int, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception7(std::function<void (int, int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        try {
            slot(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }
        catch (boost::tuple<std::string, int, int, int, int, int, int, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception8(std::function<void (int, int, int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        try {
            slot(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }
        catch (boost::tuple<std::string, int, int, int, int, int, int, int, int> const& t)
        {
            output << t;
        }
    }

    static void catch_exception9(std::function<void (int, int, int, int, int, int, int, int, int)> const& slot, int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        try {
            slot(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        }
        catch (boost::tuple<std::string, int, int, int, int, int, int, int, int, int> const& t)
        {
            output << t;
        }
    }
};

boost::test_tools::output_test_stream function_exception_fixture::output;

#endif /* ! TEST_UNIT_UTILITY_LUA_EXCEPTION_FIXTURE_HPP */
