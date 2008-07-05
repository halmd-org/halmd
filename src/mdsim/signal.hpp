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
	m_sigh[0] = signal(SIGHUP, set);
	m_sigh[1] = signal(SIGINT, set);
	m_sigh[2] = signal(SIGTERM, set);
	m_sigh[3] = signal(SIGUSR1, set);
    }

    /**
     * restore previous signal handlers
     */
    ~signal_handler()
    {
	signal(SIGHUP, m_sigh[0]);
	signal(SIGINT, m_sigh[1]);
	signal(SIGTERM, m_sigh[2]);
	signal(SIGUSR1, m_sigh[3]);
    }

    /**
     * get signal number
     */
    static int const& get()
    {
	return m_signum;
    }

    /**
     * reset signal number
     */
    void clear()
    {
	m_signum = 0;
    }

    /**
     * output signal description to stream
     */
    friend std::ostream& operator<<(std::ostream& os, signal_handler const& sig)
    {
	if (sig.m_signum == SIGHUP)
	    os << "HUP";
	else if (sig.m_signum == SIGINT)
	    os << "INT";
	else if (sig.m_signum == SIGTERM)
	    os << "TERM";
	else if (sig.m_signum == SIGUSR1)
	    os << "USR1";
	return os;
    }

private:
    /**
     * signal handler callback
     */
    static void set(int signum)
    {
	m_signum = signum;
    }

private:
    /** signal number */
    static int m_signum;
    /** previous signal handlers */
    boost::array<sighandler_t, 4> m_sigh;
};

int signal_handler::m_signum(0);

} // namespace mdsim

#endif /* ! MDSIM_SIGNAL_HPP */
