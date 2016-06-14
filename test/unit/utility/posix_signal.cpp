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

#define BOOST_TEST_MODULE posix_signal
#include <boost/test/unit_test.hpp>

#include <boost/test/parameterized_test.hpp>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

#include <halmd/utility/posix_signal.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>
#include <boost/version.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

class signal_handler
{
public:
    signal_handler(int signum) : signum_(signum), count_(0) {}

    void operator()(int signum)
    {
        BOOST_CHECK_EQUAL( signum, signum_ );
        ++count_;
    }

    size_t count() const
    {
        return count_;
    }

private:
    int signum_;
    size_t count_;
};

class test_signal_wait
{
public:
    test_signal_wait(size_t count) : count_(count) {}

    void operator()(int signum) const
    {
        BOOST_TEST_MESSAGE("wait " << posix_signal::name(signum) << ", " << count_ << " iterations");
        posix_signal sig;
        posix_signal const& sig_const(sig);
        signal_handler handler(signum);
        sig.on_signal(signum, std::ref(handler));
        for (size_t j = 0; j < count_; ++j) {
            kill(getpid(), signum);
            sig_const.wait();
        }
        BOOST_CHECK_EQUAL( handler.count(), count_ );
    }

private:
    size_t count_;
};

class test_signal_poll
{
public:
    test_signal_poll(size_t count) : count_(count) {}

    void operator()(int signum) const
    {
        BOOST_TEST_MESSAGE("poll " << posix_signal::name(signum) << ", " << count_ << " iterations");
        posix_signal sig;
        posix_signal const& sig_const(sig);
        signal_handler handler(signum);
        sig.on_signal(signum, std::ref(handler));
        for (size_t j = 0; j < count_; ++j) {
            BOOST_CHECK( !sig_const.poll() );
            kill(getpid(), signum);
            BOOST_CHECK( sig_const.poll() );
        }
        BOOST_CHECK_EQUAL( handler.count(), count_ );
    }

private:
    size_t count_;
};

BOOST_AUTO_TEST_CASE( test_alarm )
{
    BOOST_TEST_MESSAGE("wait alarm(1)");
    posix_signal sig;
    posix_signal const& sig_const(sig);
    signal_handler handler(SIGALRM);
    sig.on_signal(SIGALRM, std::ref(handler));
    timer timer;
    alarm(1);
    BOOST_CHECK( !sig_const.poll() );
    BOOST_CHECK_EQUAL( handler.count(), 0LU );
    sig_const.wait();
    BOOST_CHECK( !sig_const.poll() );
    BOOST_CHECK_CLOSE_FRACTION( timer.elapsed(), 1., 0.05 ); // 50 ms tolerance
    BOOST_CHECK_EQUAL( handler.count(), 1LU );
}

// FIXME test multiple concurrent posix_signal instances

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test;
    using namespace boost::unit_test::framework;

    vector<int> const signum = {
        SIGHUP
      , SIGINT
      , SIGALRM
      , SIGTERM
      , SIGUSR1
      , SIGUSR2
      , SIGCONT
      , SIGTSTP
      , SIGTTIN
      , SIGTTOU
    };

    for (size_t i = 1, j = 1, t; j <= 144; t = j, j += i, i = t)
    {
#if BOOST_VERSION >= 105900
        typedef boost::function<void(int)> callback_type;
#else
        typedef callback1<int> callback_type;
#endif
        master_test_suite().add( BOOST_PARAM_TEST_CASE(
            callback_type(test_signal_wait(j)), signum.begin(), signum.end()
        ) );
        master_test_suite().add( BOOST_PARAM_TEST_CASE(
            callback_type(test_signal_poll(j)), signum.begin(), signum.end()
        ) );
    }
}
