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

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/utility/lua/lua.hpp>
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
connection posix_signal::on_signal(int signum, slot_function_type const& slot)
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
    LOG_WARNING("process received signal " << name(signum));
    it->second(signum);
}

template <int signum>
static connection
wrap_on_signal(posix_signal& self, signal<void ()>::slot_function_type const& slot)
{
    return self.on_signal(signum, bind(&signal<void ()>::slot_function_type::operator(), slot));
}

static signal<void ()>::slot_function_type
wrap_wait(shared_ptr<posix_signal> self)
{
    return bind(&posix_signal::wait, self);
}

static signal<void ()>::slot_function_type
wrap_poll(shared_ptr<posix_signal> self)
{
    return bind(&posix_signal::poll, self);
}

static void
abort(shared_ptr<mdsim::clock const> clock)
{
    throw runtime_error("gracefully aborting simulation at step " + lexical_cast<string>(clock->step()));
}

static signal<void ()>::slot_function_type
wrap_abort(shared_ptr<mdsim::clock const> clock)
{
    return bind(&abort, clock);
}

void posix_signal::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("utility")
        [
            class_<posix_signal, shared_ptr<posix_signal> >("posix_signal")
                .def(constructor<>())
                .def("on_hup", &wrap_on_signal<SIGHUP>)
                .def("on_int", &wrap_on_signal<SIGINT>)
                .def("on_alrm", &wrap_on_signal<SIGALRM>)
                .def("on_term", &wrap_on_signal<SIGTERM>)
                .def("on_usr1", &wrap_on_signal<SIGUSR1>)
                .def("on_usr2", &wrap_on_signal<SIGUSR2>)
                .def("on_cont", &wrap_on_signal<SIGCONT>)
                .def("on_tstp", &wrap_on_signal<SIGTSTP>)
                .def("on_ttin", &wrap_on_signal<SIGTTIN>)
                .def("on_ttou", &wrap_on_signal<SIGTTOU>)
                .property("wait", &wrap_wait)
                .property("poll", &wrap_poll)

          , def("abort", &wrap_abort)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_utility_posix_signal(lua_State* L)
{
    posix_signal::luaopen(L);
    return 0;
}

} // namespace halmd
