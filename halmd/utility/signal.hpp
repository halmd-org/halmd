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
 * This class mimics a subset of the boost::signal interface. Many features
 * of boost::signal such as return values, connection management and object
 * tracking are not implemented. The benefit of using this minimal signal
 * class is performance close to that of boost::function calls.
 *
 * http://www.boost.org/doc/libs/release/doc/html/signals.html
 */
template <typename T>
class signal;

template <typename SlotFunction>
class signal_base
{
protected:
    typedef std::list<SlotFunction> slot_type;
    typedef boost::shared_ptr<slot_type> slot_pointer;
    typedef typename slot_type::iterator slot_iterator;
    typedef typename slot_type::const_iterator slot_const_iterator;

    slot_pointer slots_;

public:
    typedef SlotFunction slot_function_type;

    class connection
    {
    public:
        void disconnect()
        {
            slot_pointer slots = slots_.lock();
            if (slots) {
                slots->erase(iter_);
            }
        }

    private:
        friend class signal_base;

        connection(slot_pointer slots, slot_iterator iter) : slots_(slots), iter_(iter) {}

        boost::weak_ptr<slot_type> slots_;
        slot_iterator iter_;
    };

    signal_base() : slots_(new slot_type) {}

    connection connect(slot_function_type const& slot)
    {
        return connection(slots_, slots_->insert(slots_->end(), slot));
    }

    void disconnect_all_slots()
    {
        slots_->clear();
    }

    bool empty() const
    {
        return slots_->empty();
    }

    std::size_t num_slots() const
    {
        return slots_->size();
    }
};

template <typename T0>
class signal0
  : public signal_base<boost::function0<T0> >
{
public:
    T0 operator()() const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)();
        }
    }

protected:
    typedef typename signal0::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1>
class signal1
  : public signal_base<boost::function1<T0, T1> >
{
public:
    T0 operator()(T1 arg1) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1);
        }
    }

protected:
    typedef typename signal1::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2>
class signal2
  : public signal_base<boost::function2<T0, T1, T2> >
{
public:
    T0 operator()(T1 arg1, T2 arg2) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2);
        }
    }

protected:
    typedef typename signal2::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2, typename T3>
class signal3
  : public signal_base<boost::function3<T0, T1, T2, T3> >
{
public:
    T0 operator()(T1 arg1, T2 arg2, T3 arg3) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2, arg3);
        }
    }

protected:
    typedef typename signal3::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2, typename T3, typename T4>
class signal4
  : public signal_base<boost::function4<T0, T1, T2, T3, T4> >
{
public:
    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4);
        }
    }

protected:
    typedef typename signal4::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
class signal5
  : public signal_base<boost::function5<T0, T1, T2, T3, T4, T5> >
{
public:
    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5);
        }
    }

protected:
    typedef typename signal5::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
class signal6
  : public signal_base<boost::function6<T0, T1, T2, T3, T4, T5, T6> >
{
public:
    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6);
        }
    }

protected:
    typedef typename signal6::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
class signal7
  : public signal_base<boost::function7<T0, T1, T2, T3, T4, T5, T6, T7> >
{
public:
    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7);
        }
    }

protected:
    typedef typename signal7::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
class signal8
  : public signal_base<boost::function8<T0, T1, T2, T3, T4, T5, T6, T7, T8> >
{
public:
    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8);
        }
    }

protected:
    typedef typename signal8::slot_const_iterator slot_const_iterator;
};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
class signal9
  : public signal_base<boost::function9<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> >
{
public:
    T0 operator()(T1 arg1, T2 arg2, T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7, T8 arg8, T9 arg9) const
    {
        for (slot_const_iterator f = this->slots_->begin(); f != this->slots_->end(); ++f) {
            (*f)(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
        }
    }

protected:
    typedef typename signal9::slot_const_iterator slot_const_iterator;
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
