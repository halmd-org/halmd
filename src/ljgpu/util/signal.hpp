/* Signal handling
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_UTIL_SIGNAL_HPP
#define LJGPU_UTIL_SIGNAL_HPP

#include <exception>
#include <signal.h>
#include <stdlib.h>

namespace ljgpu { namespace signal
{

//
// Based on »C++ exception-handling tricks for Linux«,
// http://www.ibm.com/developerworks/linux/library/l-cppexcep.html
//

template <typename exception>
class translator
{
private:
    class singleton_translator
    {
    public:
	singleton_translator()
	{
	    struct sigaction act;
	    act.sa_handler = handler;
	    act.sa_flags = 0;
	    sigemptyset(&act.sa_mask);
	    sigaction(exception::signal(), &act, NULL);
	}

	static void handler(int)
	{
	    throw exception();
	}
    };

public:
    translator()
    {
	static singleton_translator translator;
    }
};

class HUP : public std::exception
{
public:
    static int signal() { return SIGHUP; }
    virtual char const* what() const throw() { return "SIGHUP"; }
};

class TERM : public std::exception
{
public:
    static int signal() { return SIGTERM; }
    virtual char const* what() const throw() { return "SIGTERM"; }
};

class INT : public std::exception
{
public:
    static int signal() { return SIGINT; }
    virtual char const* what() const throw() { return "SIGINT"; }
};

class ALRM : public std::exception
{
public:
    static int signal() { return SIGALRM; }
    virtual char const* what() const throw() { return "SIGALRM"; }
};

class USR1 : public std::exception
{
public:
    static int signal() { return SIGUSR1; }
    virtual char const* what() const throw() { return "SIGUSR1"; }
};

class USR2 : public std::exception
{
public:
    static int signal() { return SIGUSR2; }
    virtual char const* what() const throw() { return "SIGUSR2"; }
};

}} // namespace ljgpu::signal

#endif /* ! LJGPU_UTIL_SIGNAL_HPP */
