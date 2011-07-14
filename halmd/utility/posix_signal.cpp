/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_MDSIM_HPP
#define HALMD_MDSIM_HPP

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <halmd/utility/posix_signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {

// FIXME error code handling (with boost::system::system_error)

posix_signal::posix_signal()
{
    sigemptyset(&set_);
    pthread_sigmask(0, NULL, &oldset_);
}

posix_signal::~posix_signal()
{
    pthread_sigmask(SIG_SETMASK, &oldset_, NULL);
}

/**
 * register signal handler
 *
 * @param signum signal number (> 0)
 * @param slot signal handler function or functor
 */
posix_signal::connection_type
posix_signal::on_signal(int signum, slot_function_type const& slot)
{
    handler_map_type::iterator it;
    bool inserted;
    tie(it, inserted) = handler_.insert(make_pair(signum, handler_type()));
    if (inserted) {
        sigaddset(&set_, signum);
        pthread_sigmask(SIG_BLOCK, &set_, NULL);
    }
    return it->second.connect(slot);
}

/**
 * block process until signal is received
 */
void posix_signal::wait() const
{
    int signum = sigwaitinfo(&set_, NULL);
    this->handle(signum);
}

/**
 * poll signal queue
 *
 * @returns true if a signal was handled, false otherwise
 */
bool posix_signal::poll() const
{
    timespec timeout = { 0, 0 };
    int signum = sigtimedwait(&set_, NULL, &timeout);
    if (signum > 0) {
        this->handle(signum);
        return true;
    }
    return false;
}

/**
 * returns signal name, or signal number otherwise
 *
 * @param signum signal number (> 0)
 */
string posix_signal::name(int signum)
{
    switch (signum) {
      case SIGHUP:
        return "SIGHUP";
      case SIGINT:
        return "SIGINT";
      case SIGQUIT:
        return "SIGQUIT";
      case SIGALRM:
        return "SIGALRM";
      case SIGTERM:
        return "SIGTERM";
      case SIGUSR1:
        return "SIGUSR1";
      case SIGUSR2:
        return "SIGUSR2";
      case SIGCONT:
        return "SIGCONT";
      case SIGTSTP:
        return "SIGTSTP";
      case SIGTTIN:
        return "SIGTTIN";
      case SIGTTOU:
        return "SIGTTOU";
      default:
        return lexical_cast<string>(signum);
    }
}

/**
 * handle POSIX signal
 *
 * @param signum signal number (> 0)
 */
void posix_signal::handle(int signum) const
{
    handler_map_type::const_iterator it = handler_.find(signum);
    if (it == handler_.end()) {
        throw std::logic_error("blocked unregistered signal " + name(signum));
    }
    it->second(signum);
}

} // namespace halmd

#endif /* ! HALMD_MDSIM_HPP */
