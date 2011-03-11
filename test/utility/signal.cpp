/*
 * Copyright © 2011  Peter Colberg
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

#define BOOST_TEST_MODULE signal
#include <boost/test/unit_test.hpp>

#include <boost/bind.hpp>

#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

class signal_counter_base
{
protected:
    size_t count_;

public:
    signal_counter_base() : count_(0) {}

    size_t count() const
    {
        return count_;
    }
};

class signal_counter0
  : public signal_counter_base
{
public:
    void operator()()
    {
        count_ += 1LU;
    }
};

class signal_counter1
  : public signal_counter_base
{
public:
    void operator()(int arg1)
    {
        count_ += 1LU;
        count_ += arg1;
    }
};

class signal_counter2
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
    }
};

class signal_counter3
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2, int arg3)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
        count_ += arg3;
    }
};

class signal_counter4
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2, int arg3, int arg4)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
        count_ += arg3;
        count_ += arg4;
    }
};

class signal_counter5
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2, int arg3, int arg4, int arg5)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
        count_ += arg3;
        count_ += arg4;
        count_ += arg5;
    }
};

class signal_counter6
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
        count_ += arg3;
        count_ += arg4;
        count_ += arg5;
        count_ += arg6;
    }
};

class signal_counter7
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
        count_ += arg3;
        count_ += arg4;
        count_ += arg5;
        count_ += arg6;
        count_ += arg7;
    }
};

class signal_counter8
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
        count_ += arg3;
        count_ += arg4;
        count_ += arg5;
        count_ += arg6;
        count_ += arg7;
        count_ += arg8;
    }
};

class signal_counter9
  : public signal_counter_base
{
public:
    void operator()(int arg1, int arg2, int arg3, int arg4, int arg5, int arg6, int arg7, int arg8, int arg9)
    {
        count_ += 1LU;
        count_ += arg1;
        count_ += arg2;
        count_ += arg3;
        count_ += arg4;
        count_ += arg5;
        count_ += arg6;
        count_ += arg7;
        count_ += arg8;
        count_ += arg9;
    }
};

template <typename signal_type>
size_t result(size_t calls, size_t slots)
{
    size_t count = signal_type::slot_function_type::arity + 1;
    return (count * (count + 1) / 2) * calls * slots;
}

BOOST_AUTO_TEST_CASE( halmd_signal0 )
{
    typedef signal<void ()> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter0 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig();
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig();
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig();
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal1 )
{
    typedef signal<void (int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter1 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal2 )
{
    typedef signal<void (int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter2 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal3 )
{
    typedef signal<void (int, int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter3 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3, 4);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3, 4);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3, 4);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal4 )
{
    typedef signal<void (int, int, int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter4 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), 0LU );
    BOOST_CHECK_EQUAL( counter2.count(), 0LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3, 4, 5);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3, 4, 5);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3, 4, 5);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal5 )
{
    typedef signal<void (int, int, int, int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter5 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3, 4, 5, 6);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3, 4, 5, 6);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3, 4, 5, 6);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal6 )
{
    typedef signal<void (int, int, int, int, int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter6 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3, 4, 5, 6, 7);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3, 4, 5, 6, 7);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3, 4, 5, 6, 7);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal7 )
{
    typedef signal<void (int, int, int, int, int, int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter7 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), 0LU );
    BOOST_CHECK_EQUAL( counter2.count(), 0LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3, 4, 5, 6, 7, 8);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3, 4, 5, 6, 7, 8);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3, 4, 5, 6, 7, 8);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal8 )
{
    typedef signal<void (int, int, int, int, int, int, int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter8 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3, 4, 5, 6, 7, 8, 9);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3, 4, 5, 6, 7, 8, 9);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3, 4, 5, 6, 7, 8, 9);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}

BOOST_AUTO_TEST_CASE( halmd_signal9 )
{
    typedef signal<void (int, int, int, int, int, int, int, int, int)> signal_type;

    signal_type sig;
    signal_type const& immutable_sig(sig);
    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    signal_counter9 counter1, counter2;
    sig.connect(ref(counter1));
    sig.connect(ref(counter1));
    signal_type::connection c = sig.connect(ref(counter2));
    sig.connect(ref(counter1));

    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 4LU );
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(0, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(0, 1) );

    immutable_sig(2, 3, 4, 5, 6, 7, 8, 9, 10);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(1, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect(c);
    BOOST_CHECK( !immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 3LU );

    immutable_sig(2, 3, 4, 5, 6, 7, 8, 9, 10);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );

    sig.disconnect_all_slots();

    BOOST_CHECK( immutable_sig.empty() );
    BOOST_CHECK_EQUAL( immutable_sig.num_slots(), 0LU );

    immutable_sig(2, 3, 4, 5, 6, 7, 8, 9, 10);
    BOOST_CHECK_EQUAL( counter1.count(), result<signal_type>(2, 3) );
    BOOST_CHECK_EQUAL( counter2.count(), result<signal_type>(1, 1) );
}
