/* Signal handling
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_SIGNAL_HPP
#define MDSIM_SIGNAL_HPP

#include <boost/array.hpp>
#include <signal.h>


namespace mdsim
{

/**
 * Signal handling
 */
class signal_handler
{
public:
    /**
     * set signal handlers
     */
    signal_handler()
    {
	handler[0] = ::signal(SIGALRM, callback);
	handler[1] = ::signal(SIGHUP, callback);
	handler[2] = ::signal(SIGINT, callback);
	handler[3] = ::signal(SIGTERM, callback);
	handler[4] = ::signal(SIGUSR1, callback);
	handler[5] = ::signal(SIGUSR2, callback);
    }

    /**
     * restore previous signal handlers
     */
    ~signal_handler()
    {
	::signal(SIGALRM, handler[0]);
	::signal(SIGHUP, handler[1]);
	::signal(SIGINT, handler[2]);
	::signal(SIGTERM, handler[3]);
	::signal(SIGUSR1, handler[4]);
	::signal(SIGUSR2, handler[5]);
    }

    /**
     * get signal number
     */
    static int const& get()
    {
	return signum;
    }

    /**
     * set signal number to zero
     */
    void clear()
    {
	signum = 0;
    }

private:
    /**
     * signal handler callback
     */
    static void callback(int signum)
    {
	signal_handler::signum = signum;
    }

private:
    /** signal number */
    static int signum;
    /** previous signal handlers */
    boost::array<sighandler_t, 6> handler;
};

int signal_handler::signum(0);

} // namespace mdsim

#endif /* ! MDSIM_SIGNAL_HPP */
