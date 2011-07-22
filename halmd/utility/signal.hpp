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

#ifndef HALMD_UTILITY_SIGNAL_HPP
#define HALMD_UTILITY_SIGNAL_HPP

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <list>

namespace halmd {

/**
 * Signal implements callbacks with multiple targets.
 *
 * This class mimics a subset of the boost::signal interface. Some advanced
 * features of boost::signal such as return values and object tracking are
 * not implemented. The benefit of using this minimal signal class is
 * performance close to that of boost::function calls.
 *
 * http://www.boost.org/doc/libs/release/doc/html/signals.html
 */
template <typename T>
class signal;

/**
 * Slot to slots container connection
 *
 * The connection is a slot itself, so it may be connected to
 * a signal<void ()>, to disconnect this slot when invoked.
 */
typedef boost::function0<void> connection;

template <typename T>
class slots
{
private:
    /**
     * As storage for the slots, we use a linked list instead of a vector,
     * which guarantees that an iterator to an inserted slot remains valid
     * as long as the slot is not removed from the list, i.e. the iterator
     * remains valid after insertion or removal of other slots.
     */
    typedef std::list<T> slots_type;
    /**
     * The list of slots is held with a shared pointer, which allows
     * tracking the slots with a weak pointer in connection objects.
     */
    typedef boost::shared_ptr<slots_type> slots_pointer;
    /**
     * We do not provide write access to the list container,
     * to ensure that connection iterators remain valid,
     * therefore declare iterator private.
     */
    typedef typename slots_type::iterator iterator;

    /**
     * Slot-to-signal connection
     *
     * This class manages the connection of a single slot to a signal
     * *independent* of the lifetime of the signal. We employ a weak
     * pointer to observe the shared pointer holding the slots in the
     * signal object.
     * A weak pointer is locked to retrieve a shared pointer to the
     * referenced object. If the weak pointer does not observe an
     * object, lock() will return an empty shared pointer instead.
     * This is the case if (1) we erase() the slot using the stored
     * iterator and reset() the weak pointer, (2) all slots are
     * disconnected from the signal with disconnect_all_slots()
     * by reallocating the slots list, or (3) the signal object
     * along with the slots list is deconstructed. Thus the state of
     * the weak pointer reflects the validity of the connection object.
     */
    class connection_
    {
    public:
        /**
         * disconnect slot from slots container
         *
         * This method may be invoked multiple times. Repeated calls to
         * disconnect() will be silently ignored, similar to the behaviour
         * of boost::signals::connection::disconnect(). In particular, this
         * method may be called even if the slots container no longer exists,
         * or if all slots have been removed with disconnect_all_slots().
         */
        void operator()()
        {
            slots_pointer slots = slots_.lock();
            if (slots) {
                slots->erase(iter_);
                slots_.reset();
            }
        }

        connection_(slots_pointer slots, iterator iter) : slots_(slots), iter_(iter) {}

    private:
        boost::weak_ptr<slots_type> slots_;
        iterator iter_;
    };

    slots_pointer slots_;

public:
    /**
     * We provide read access to the list container with slots::begin()
     * and slots::end(), therefore declare const_iterator public.
     */
    typedef typename slots_type::const_iterator const_iterator;

    slots() : slots_(new slots_type) {}

    /**
     * connect slot to signal
     */
    connection connect(T const& slot)
    {
        return connection_(slots_, slots_->insert(slots_->end(), slot));
    }

    /**
     * disconnect all slots
     *
     * Instead of clearing the list of slots, we reallocate the list,
     * which breaks the link between weak pointers in connection objects
     * and the shared_ptr in signal holding the slots, and enables
     * connection objects to ignore calls to connection::disconnect().
     */
    void disconnect_all_slots()
    {
        slots_.reset(new slots_type);
    }

    /**
     * returns true if list of slots is empty, false otherwise
     */
    bool empty() const
    {
        return slots_->empty();
    }

    /**
     * returns number of connected slots
     */
    std::size_t num_slots() const
    {
        return slots_->size();
    }

    /**
     * returns const iterator to first slot
     */
    const_iterator begin() const
    {
        return slots_->begin();
    }

    /**
     * returns const iterator one past last slot
     */
    const_iterator end() const
    {
        return slots_->end();
    }
};

template <typename T0>
class signal0
  : public slots<boost::function0<T0> >
{
public:
    typedef boost::function0<T0> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()() const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)();
        }
    }
};

template <typename T0, typename T1>
class signal1
  : public slots<boost::function1<T0, T1> >
{
public:
    typedef boost::function1<T0, T1> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1);
        }
    }
};

template <typename T0, typename T1, typename T2>
class signal2
  : public slots<boost::function2<T0, T1, T2> >
{
public:
    typedef boost::function2<T0, T1, T2> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2);
        }
    }
};

template <typename T0, typename T1, typename T2, typename T3>
class signal3
  : public slots<boost::function3<T0, T1, T2, T3> >
{
public:
    typedef boost::function3<T0, T1, T2, T3> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2, T3 arg3) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3);
        }
    }
};

template <typename T0, typename T1, typename T2, typename T3, typename T4>
class signal4
  : public slots<boost::function4<T0, T1, T2, T3, T4> >
{
public:
    typedef boost::function4<T0, T1, T2, T3, T4> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4);
        }
    }
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
class signal5
  : public slots<boost::function5<T0, T1, T2, T3, T4, T5> >
{
public:
    typedef boost::function5<T0, T1, T2, T3, T4, T5> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5);
        }
    }
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
class signal6
  : public slots<boost::function6<T0, T1, T2, T3, T4, T5, T6> >
{
public:
    typedef boost::function6<T0, T1, T2, T3, T4, T5, T6> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6);
        }
    }
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
class signal7
  : public slots<boost::function7<T0, T1, T2, T3, T4, T5, T6, T7> >
{
public:
    typedef boost::function7<T0, T1, T2, T3, T4, T5, T6, T7> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }
    }
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
class signal8
  : public slots<boost::function8<T0, T1, T2, T3, T4, T5, T6, T7, T8> >
{
public:
    typedef boost::function8<T0, T1, T2, T3, T4, T5, T6, T7, T8> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }
    }
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
class signal9
  : public slots<boost::function9<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> >
{
public:
    typedef boost::function9<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8, T9 arg9) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        }
    }
};

template <>
class signal<void ()>
  : public signal0<void> {};

template <typename T1>
class signal<void (T1)>
  : public signal1<void, T1> {};

template <typename T1, typename T2>
class signal<void (T1, T2)>
  : public signal2<void, T1, T2> {};

template <typename T1, typename T2, typename T3>
class signal<void (T1, T2, T3)>
  : public signal3<void, T1, T2, T3> {};

template <typename T1, typename T2, typename T3, typename T4>
class signal<void (T1, T2, T3, T4)>
  : public signal4<void, T1, T2, T3, T4> {};

template <typename T1, typename T2, typename T3, typename T4, typename T5>
class signal<void (T1, T2, T3, T4, T5)>
  : public signal5<void, T1, T2, T3, T4, T5> {};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
class signal<void (T1, T2, T3, T4, T5, T6)>
  : public signal6<void, T1, T2, T3, T4, T5, T6> {};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
class signal<void (T1, T2, T3, T4, T5, T6, T7)>
  : public signal7<void, T1, T2, T3, T4, T5, T6, T7> {};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
class signal<void (T1, T2, T3, T4, T5, T6, T7, T8)>
  : public signal8<void, T1, T2, T3, T4, T5, T6, T7, T8> {};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
class signal<void (T1, T2, T3, T4, T5, T6, T7, T8, T9)>
  : public signal9<void, T1, T2, T3, T4, T5, T6, T7, T8, T9> {};

} // namespace halmd

#endif /* ! HALMD_UTILITY_SIGNAL_HPP */
