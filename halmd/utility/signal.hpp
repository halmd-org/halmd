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

// forward declaration
template <typename T>
class slots;

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
class connection
{
private:
    typedef boost::shared_ptr<void> slots_pointer;
    typedef boost::function<void (slots_pointer)> disconnector_type;

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
    void disconnect()
    {
        slots_pointer slots = slots_.lock();
        if (slots) {
            disconnect_(slots);
            slots_.reset();
        }
    }

    /**
     * returns true if slot is connected to slots container, false otherwise
     */
    bool connected() const
    {
        return slots_.lock();
    }

private:
    /**
     * Only a slots object should be able to construct a connection.
     */
    template <typename T>
    friend class slots;

    connection(slots_pointer slots, disconnector_type c) : slots_(slots), disconnect_(c) {}

    boost::weak_ptr<void> slots_;
    disconnector_type disconnect_;
};

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
    typedef typename slots_type::iterator slots_iterator;

    /**
     * This class disconnects the stored slot from the slots container.
     * It is wrapped as a type-independent function object for use in
     * a connection object, which triggers the disconnect only if the
     * slot and the slots container still exist.
     */
    class disconnector
    {
    public:
        void operator()(boost::shared_ptr<void> slots)
        {
            slots_pointer slots_ = boost::static_pointer_cast<slots_type>(slots);
            slots_->erase(iter_);
        }

        disconnector(slots_iterator iter) : iter_(iter) {}

    private:
        slots_iterator iter_;
    };

    slots_pointer slots_;

public:
    /**
     * We provide read access to the list container with slots::begin()
     * and slots::end(), therefore declare const_iterator public.
     */
    typedef typename slots_type::const_iterator const_iterator;
    /**
     * Alias const_iterator as iterator for use with BOOST_FOREACH.
     *
     * http://www.boost.org/doc/libs/release/doc/html/foreach/extensibility.html
     * http://www.boost.org/doc/libs/release/libs/range/doc/html/range/reference/extending/method_1.html
     */
    typedef const_iterator iterator;

    slots() : slots_(new slots_type) {}

    /**
     * connect slot to signal
     */
    connection connect(T const& slot)
    {
        slots_iterator iter = slots_->insert(slots_->end(), slot);
        return connection(slots_, disconnector(iter));
    }

    /**
     * disconnect all slots
     *
     * Instead of clearing the list of slots, we reallocate the list,
     * which breaks the link between weak pointers in connection objects
     * and the boost::shared_ptr in signal holding the slots, and enables
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
template <typename Function>
class signal;

#ifndef HALMD_NO_CXX11

template <typename... Args>
class signal<void (Args...)>
  : public slots<boost::function<void (Args...)> >
{
public:
    typedef boost::function<void (Args...)> slot_function_type;

    void operator()(Args... args) const
    {
        for (slot_function_type const& f : *this) {
            f(args...);
        }
    }
};

#else /* HALMD_NO_CXX11 */

template <>
class signal<void ()>
  : public slots<boost::function<void ()> >
{
public:
    typedef boost::function<void ()> slot_function_type;
    typedef slots<slot_function_type>::const_iterator const_iterator;

    void operator()() const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)();
        }
    }
};

template <typename T1>
class signal<void (T1)>
  : public slots<boost::function<void (T1)> >
{
public:
    typedef boost::function<void (T1)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1);
        }
    }
};

template <typename T1, typename T2>
class signal<void (T1, T2)>
  : public slots<boost::function<void (T1, T2)> >
{
public:
    typedef boost::function<void (T1, T2)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2);
        }
    }
};

template <typename T1, typename T2, typename T3>
class signal<void (T1, T2, T3)>
  : public slots<boost::function<void (T1, T2, T3)> >
{
public:
    typedef boost::function<void (T1, T2, T3)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2, T3 arg3) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3);
        }
    }
};

template <typename T1, typename T2, typename T3, typename T4>
class signal<void (T1, T2, T3, T4)>
  : public slots<boost::function<void (T1, T2, T3, T4)> >
{
public:
    typedef boost::function<void (T1, T2, T3, T4)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4);
        }
    }
};

template <typename T1, typename T2, typename T3, typename T4, typename T5>
class signal<void (T1, T2, T3, T4, T5)>
  : public slots<boost::function<void (T1, T2, T3, T4, T5)> >
{
public:
    typedef boost::function<void (T1, T2, T3, T4, T5)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5);
        }
    }
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
class signal<void (T1, T2, T3, T4, T5, T6)>
  : public slots<boost::function<void (T1, T2, T3, T4, T5, T6)> >
{
public:
    typedef boost::function<void (T1, T2, T3, T4, T5, T6)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6);
        }
    }
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
class signal<void (T1, T2, T3, T4, T5, T6, T7)>
  : public slots<boost::function<void (T1, T2, T3, T4, T5, T6, T7)> >
{
public:
    typedef boost::function<void (T1, T2, T3, T4, T5, T6, T7)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }
    }
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
class signal<void (T1, T2, T3, T4, T5, T6, T7, T8)>
  : public slots<boost::function<void (T1, T2, T3, T4, T5, T6, T7, T8)> >
{
public:
    typedef boost::function<void (T1, T2, T3, T4, T5, T6, T7, T8)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }
    }
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
class signal<void (T1, T2, T3, T4, T5, T6, T7, T8, T9)>
  : public slots<boost::function<void (T1, T2, T3, T4, T5, T6, T7, T8, T9)> >
{
public:
    typedef boost::function<void (T1, T2, T3, T4, T5, T6, T7, T8, T9)> slot_function_type;
    typedef typename slots<slot_function_type>::const_iterator const_iterator;

    void operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8, T9 arg9) const
    {
        for (const_iterator f = this->begin(); f != this->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        }
    }
};

#endif /* HALMD_NO_CXX11 */

} // namespace halmd

#endif /* ! HALMD_UTILITY_SIGNAL_HPP */
