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

#include <algorithm>
#include <signal.h>
#include <stdlib.h>

namespace ljgpu { namespace signal
{

/**
 * synchronous signal handler
 */
class handler
{
public:
    /**
     * block all signals in this process
     */
    handler() {
        sigset_t set;
        sigfillset(&set);
        sigprocmask(SIG_BLOCK, &set, NULL);
    }

    /**
     * unblock all signals in this process
     */
    ~handler() {
        sigset_t set;
        sigfillset(&set);
        sigprocmask(SIG_UNBLOCK, &set, NULL);
    }
};

extern int signal;

/**
 * poll signal queue
 */
inline int poll()
{
    timespec tv;
    tv.tv_sec = 0;
    tv.tv_nsec = 0;
    sigset_t set;
    sigfillset(&set);
    signal = std::max(0, sigtimedwait(&set, NULL, &tv));
    return signal;
}

/**
 * block process until signal is received
 */
inline int wait()
{
    sigset_t set;
    sigfillset(&set);
    signal = std::max(0, sigwaitinfo(&set, NULL));
    return signal;
}

}} // namespace ljgpu::signal

#endif /* ! LJGPU_UTIL_SIGNAL_HPP */
