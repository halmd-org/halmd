/*
 * Copyright © 2008-2012  Peter Colberg
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
#include <boost/system/error_code.hpp>
#include <boost/system/system_error.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/posix_signal.hpp>

using namespace std;

namespace halmd {

/**
 * Block all signals
 *
 * @returns signal set with blocked signals
 */
sigset_t posix_signal::block_signals()
{
    namespace bs = boost::system;

    sigset_t set;
    sigfillset(&set);
    bs::error_code ec(pthread_sigmask(SIG_SETMASK, &set, NULL), bs::get_posix_category());
    if (ec != bs::errc::success) {
        throw bs::system_error(ec);
    }
    return set;
}

/**
 * Block all signals on program startup
 *
 * The CUDA library interfers with signal blocking.
 *
 * We can only speculate as to what is happening in libcuda.so:
 *
 *    # strings /usr/lib64/libcuda.so | grep ^sig
 *    sigemptyset
 *    sigaddset
 *    sigprocmask
 *
 * If we block signals before creating a CUDA context in the device module,
 * signal handling works as expected. If we block signals after creating a
 * CUDA context, signal handling fails, e.g. on TERM the process aborts.
 *
 * As a work-around, we block *all* signals using static initialisation,
 * and register the handlers later, upon construction of posix_signal.
 *
 * The static variable posix_signal::set_ is used by class members, which
 * prevents the compiler from optimizing away the call to block_signals().
 */
sigset_t posix_signal::set_ = posix_signal::block_signals();

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
        return "HUP";
      case SIGINT:
        return "INT";
      case SIGALRM:
        return "ALRM";
      case SIGTERM:
        return "TERM";
      case SIGUSR1:
        return "USR1";
      case SIGUSR2:
        return "USR2";
      case SIGCONT:
        return "CONT";
      case SIGTSTP:
        return "TSTP";
      case SIGTTIN:
        return "TTIN";
      case SIGTTOU:
        return "TTOU";
      default:
        return boost::lexical_cast<string>(signum);
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
    if (it != handler_.end()) {
        LOG_WARNING("process received signal " << name(signum));
        it->second(signum);
    }
}

template <int signum>
static connection
wrap_on_signal(posix_signal& self, function<void ()> const& slot)
{
    return self.on_signal(signum, [=](int) {
        slot();
    });
}

static function<void ()>
wrap_wait(shared_ptr<posix_signal> self)
{
    return [=]() {
        self->wait();
    };
}

static function<void ()>
wrap_poll(shared_ptr<posix_signal> self)
{
    return [=]() {
        self->poll();
    };
}

void posix_signal::luaopen(lua_State* L)
{
    using namespace luaponte;
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
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_utility_posix_signal(lua_State* L)
{
    posix_signal::luaopen(L);
    return 0;
}

} // namespace halmd
