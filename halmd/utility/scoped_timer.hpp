/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_UTILITY_SCOPED_TIMER_HPP
#define HALMD_UTILITY_SCOPED_TIMER_HPP

#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace halmd {

/**
 * Scoped timer
 */
template <typename Timer>
class scoped_timer
{
public:
    /**
     * Start timer.
     */
    template <typename UnaryFunctor>
    scoped_timer(UnaryFunctor& f);

    /**
     * Invoke functor with elapsed time.
     */
    ~scoped_timer();

private:
    boost::function<void ()> elapsed_;

    template <typename UnaryFunctor>
    static void elapsed(UnaryFunctor& f, Timer const& t);
};

template <typename Timer>
template <typename UnaryFunctor>
inline scoped_timer<Timer>::scoped_timer(UnaryFunctor& f)
  : elapsed_(
        boost::bind(
            &scoped_timer<Timer>::template elapsed<UnaryFunctor>
          , boost::ref(f)
          , Timer() // start timer, as late as possible
        )
    ) {}

template <typename Timer>
inline scoped_timer<Timer>::~scoped_timer()
{
    elapsed_();
}

template <typename Timer>
template <typename UnaryFunctor>
inline void scoped_timer<Timer>::elapsed(UnaryFunctor& f, Timer const& t)
{
    f(t.elapsed());
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_SCOPED_TIMER_HPP */
