/*
 * Copyright © 2011  Peter Colberg
 * Copyright © 2020 Jaslo Ziska
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

#define BOOST_TEST_MODULE posix_signal
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <boost/test/parameterized_test.hpp>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

#include <halmd/utility/posix_signal.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
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
using namespace boost::unit_test;

size_t const DATA_ARRAY_COUNT[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144};
int const DATA_ARRAY_SIGNUM[] = {
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
auto dataset = data::make(DATA_ARRAY_COUNT) * data::make(DATA_ARRAY_SIGNUM);

BOOST_DATA_TEST_CASE( test_signal_wait, dataset, count, signum) {
    BOOST_TEST_MESSAGE("wait " << posix_signal::name(signum) << ", " << count << " iterations");
    posix_signal sig;
    posix_signal const& sig_const(sig);
    signal_handler handler(signum);
    sig.on_signal(signum, std::ref(handler));
    for (size_t j = 0; j < count; ++j) {
        kill(getpid(), signum);
        sig_const.wait();
    }
    BOOST_CHECK_EQUAL( handler.count(), count );
}

BOOST_DATA_TEST_CASE( test_signal_poll, dataset, count, signum) {
    BOOST_TEST_MESSAGE("poll " << posix_signal::name(signum) << ", " << count << " iterations");
    posix_signal sig;
    posix_signal const& sig_const(sig);
    signal_handler handler(signum);
    sig.on_signal(signum, std::ref(handler));
    for (size_t j = 0; j < count; ++j) {
        BOOST_CHECK( !sig_const.poll() );
        kill(getpid(), signum);
        BOOST_CHECK( sig_const.poll() );
    }
    BOOST_CHECK_EQUAL( handler.count(), count );
}
