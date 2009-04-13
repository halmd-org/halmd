/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   slim_string.cpp
 * \author Andrey Semashev
 * \date   19.11.2007
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <stdexcept>
#include <algorithm>
#include <functional>
#include <boost/log/detail/throw_exception.hpp>
#include <boost/log/utility/slim_string.hpp>
#ifndef BOOST_LOG_NO_THREADS
#include <boost/detail/atomic_count.hpp>
#endif // BOOST_LOG_NO_THREADS
#include "alignment_gap_between.hpp"

#ifdef _MSC_VER
#pragma warning(push)
// 'this' : used in base member initializer list
#pragma warning(disable: 4355)
#endif // _MSC_VER

namespace boost {

namespace BOOST_LOG_NAMESPACE {

//! Slim string implementation class
template< typename CharT, typename TraitsT >
struct basic_slim_string< CharT, TraitsT >::implementation
{
private:
    //! Auxiliary function object to implement character comparison through char_traits
    struct eq_traits :
        public std::binary_function< char_type, char_type, bool >
    {
        bool operator()(char_type left, char_type right) const
        {
            return traits_type::eq(left, right);
        }
    };
    //! Auxiliary function object to implement character lookup
    struct eq_char_bound :
        public std::unary_function< char_type, bool >
    {
    private:
        char_type m_val;

    public:
        eq_char_bound(char_type c) : m_val(c) {}
        bool operator()(char_type x) const
        {
            return traits_type::eq(x, m_val);
        }
    };

private:
    //! Reference counter
#ifndef BOOST_LOG_NO_THREADS
    detail::atomic_count m_RefCount;
#else
    unsigned long m_RefCount;
#endif // BOOST_LOG_NO_THREADS

    //! Pointer to the first character
    pointer m_pBegin;
    //! Pointer to the last character
    pointer m_pEnd;

public:
    //! Constructor
    implementation(pointer b, const_pointer str, size_type n) : m_RefCount(1), m_pBegin(b), m_pEnd(b + n)
    {
        if (str)
            traits_type::copy(b, str, n);
        b[n] = 0;
    }

    //! Returns begin iterator
    pointer begin() const { return m_pBegin; }
    //! Returns end iterator
    pointer end() const { return m_pEnd; }
    //! Returns string length
    size_type size() const { return (m_pEnd - m_pBegin); }

    //  Lookup methods
    size_type find(const_pointer that, size_type len, size_type pos) const
    {
        const size_type size = this->size();
        if (pos < size && len < (size - pos))
        {
            const_iterator p =
                std::search(m_pBegin + pos, m_pEnd, that, that + len, eq_traits());
            return (p != m_pEnd ? (size_type)(p - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type find(char_type c, size_type pos) const
    {
        const size_type size = this->size();
        if (pos < size)
        {
            register const_pointer p = traits_type::find(m_pBegin + pos, size - pos, c);
            return (p != NULL ? (size_type)(p - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type rfind(const_pointer that, size_type len, size_type pos) const
    {
        const size_type size = this->size();
        if (size >= len)
        {
            register const_pointer pBegin = m_pBegin + (std::min)(pos, size - len);
            while (pBegin >= m_pBegin)
            {
                if (0 == traits_type::compare(pBegin, that, len))
                    return (size_type)(pBegin - m_pBegin);
                --pBegin;
            }
        }
        return (size_type)npos;
    }
    size_type rfind(value_type c, size_type pos) const
    {
        const size_type size = this->size();
        if (size > 0)
        {
            register const_pointer pBegin = m_pBegin + (std::min)(pos, size - 1) + 1;
            const_reverse_iterator rend(m_pBegin);
            const_reverse_iterator r = std::find_if(const_reverse_iterator(pBegin), rend, eq_char_bound(c));
            return (r != rend ? (size_type)((r.base() - 1) - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type find_first_of(const_pointer collection, size_type collection_size, size_type pos) const
    {
        if (pos < size())
        {
            register const_pointer const p =
                _find_first_of(m_pBegin + pos, m_pEnd, collection, collection_size);
            return ((p != m_pEnd) ? (size_type)(p - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type find_last_of(const_pointer collection, size_type collection_size, size_type pos) const
    {
        const size_type size = this->size();
        if (size > 0)
        {
            register const_pointer pBegin = m_pBegin + (std::min)(pos, size - 1) + 1;
            const_reverse_iterator rend(m_pBegin);
            const_reverse_iterator r = _find_first_of(const_reverse_iterator(pBegin), rend, collection, collection_size);
            return ((r != rend) ? (size_type)((r.base() - 1) - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type find_first_not_of(const_pointer collection, size_type collection_size, size_type pos) const
    {
        if (pos < size())
        {
            register const_pointer const p =
                _find_first_not_of(m_pBegin + pos, m_pEnd, collection, collection_size);
            return ((p != m_pEnd) ? (size_type)(p - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type find_last_not_of(const_pointer collection, size_type collection_size, size_type pos) const
    {
        const size_type size = this->size();
        if (size > 0)
        {
            register const_pointer pBegin = m_pBegin + (std::min)(pos, size - 1) + 1;
            const_reverse_iterator rend(m_pBegin);
            const_reverse_iterator r = _find_first_not_of(const_reverse_iterator(pBegin), rend, collection, collection_size);
            return ((r != rend) ? (size_type)((r.base() - 1) - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type find_first_not_of(value_type c, size_type pos) const
    {
        if (pos < size())
        {
            register const_pointer p = m_pBegin + pos;
            while (p < m_pEnd && traits_type::eq(c, *p))
                ++p;
            return (p != m_pEnd ? (size_type)(p - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }
    size_type find_last_not_of(value_type c, size_type pos) const
    {
        const size_type size = this->size();
        if (size > 0)
        {
            register const_pointer p = m_pBegin + (std::min)(pos, size - 1);
            while (p >= m_pBegin && traits_type::eq(c, *p))
                --p;
            return (p >= m_pBegin ? (size_type)(p - m_pBegin) : (size_type)npos);
        }
        else
            return (size_type)npos;
    }

    //! Comparison implementation
    int compare(size_type pos1, size_type n1, const_pointer that, size_type len) const
    {
        const size_type size = this->size();
        if (pos1 <= size)
        {
            const size_type part_size = (std::min)(n1, (size_type)size - pos1);
            register const int res = traits_type::compare(m_pBegin + pos1, that, (std::min)(part_size, len));
            if (res == 0)
                return static_cast< int >(part_size - len);
            else
                return res;
        }
        else
        {
            boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::compare: the position is out of range"));
        }
    }

    //! Increments the reference counter
    void add_ref() { ++m_RefCount; }
    //! Decrements the reference counter
    unsigned long release() { return static_cast< unsigned long >(--m_RefCount); }

    //! Allocates the needed memory to accomodate n characters and constructs implementation object with m_RefCount == 1.
    static implementation* allocate_and_add_ref(internal_allocator_type* pAlloc, const_pointer str, size_type n)
    {
        enum { alignment_gap = log::aux::alignment_gap_between< implementation, char_type >::value };
        typename internal_allocator_type::pointer p =
            pAlloc->allocate(sizeof(implementation) + ((n + 1) * sizeof(char_type)) + alignment_gap);
        new (p) implementation(reinterpret_cast< char_type* >(p + sizeof(implementation) + alignment_gap), str, n);
        return reinterpret_cast< implementation* >(p);
    }
    //! Releases the implementation object and if its ref. counter drops to 0, destroys it and deallocates memory.
    static void release_and_dispose(implementation* pImpl, internal_allocator_type* pAlloc)
    {
        if (pImpl->release() == 0)
        {
            enum { alignment_gap = log::aux::alignment_gap_between< implementation, char_type >::value };
            const size_type size =
                sizeof(implementation) + ((pImpl->size() + 1) * sizeof(char_type)) + alignment_gap;
            pImpl->~implementation();
            pAlloc->deallocate(
                reinterpret_cast< typename internal_allocator_type::pointer >(pImpl),
                size);
        }
    }

private:
    template< typename IteratorT >
    static BOOST_LOG_FORCEINLINE IteratorT _find_first_of(
        IteratorT first, IteratorT last, const_pointer collection, size_type collection_size)
    {
        for (; first != last; ++first)
        {
            for (register size_type i = collection_size; i > 0; --i)
            {
                if (traits_type::eq(*first, collection[i - 1]))
                    return first;
            }
        }
        return last;
    }
    template< typename IteratorT >
    static BOOST_LOG_FORCEINLINE IteratorT _find_first_not_of(
        IteratorT first, IteratorT last, const_pointer collection, size_type collection_size)
    {
        for (; first != last; ++first)
        {
            register size_type i = collection_size;
            for (; i > 0; --i)
            {
                if (traits_type::eq(*first, collection[i - 1]))
                    break;
            }
            if (i == 0)
                return first;
        }
        return last;
    }
};


//! Default constructor
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string() : m_pImpl(implementation::allocate_and_add_ref(this, NULL, 0)) {}
//! Copy constructor
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string(basic_slim_string const& that) : m_pImpl(that.m_pImpl)
{
    m_pImpl->add_ref();
}

//  Standard constructors
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string(string_type const& that)
    : m_pImpl(implementation::allocate_and_add_ref(this, that.data(), that.size()))
{
}
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string(string_type const& s, size_type pos, size_type n)
{
    const size_type that_size = s.size();
    if (pos <= that_size)
    {
        m_pImpl = implementation::allocate_and_add_ref(this, s.data() + pos, (std::min)(n, that_size - pos));
    }
    else
    {
        boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::basic_slim_string: the position is out of range"));
    }
}
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string(basic_slim_string const& that, size_type pos, size_type n)
{
    const size_type that_size = that.m_pImpl->size();
    if (pos == 0 && n >= that_size)
    {
        m_pImpl = that.m_pImpl;
        m_pImpl->add_ref();
    }
    else
    {
        if (pos <= that_size)
        {
            m_pImpl = implementation::allocate_and_add_ref(this, that.data() + pos, (std::min)(n, that_size - pos));
        }
        else
        {
            boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::basic_slim_string: the position is out of range"));
        }
    }
}
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string(const_pointer s)
    : m_pImpl(implementation::allocate_and_add_ref(this, s, traits_type::length(s))) {}
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string(const_pointer s, size_type n)
    : m_pImpl(implementation::allocate_and_add_ref(this, s, n)) {}
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::basic_slim_string(size_type n, char_type c)
    : m_pImpl(implementation::allocate_and_add_ref(this, NULL, n))
{
    traits_type::assign(m_pImpl->begin(), n, c);
}

//! Destructor
template< typename CharT, typename TraitsT >
basic_slim_string< CharT, TraitsT >::~basic_slim_string()
{
    implementation::release_and_dispose(m_pImpl, this);
}

//! Indexing
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::const_reference
#else
CharT const&
#endif
basic_slim_string< CharT, TraitsT >::operator[] (size_type n) const { return *(m_pImpl->begin() + n); }

//  Accessors
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::const_reference
#else
CharT const&
#endif
basic_slim_string< CharT, TraitsT >::at(size_type n) const
{
    if (n >= m_pImpl->size())
        boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::at: character index is out of range"));

    return *(m_pImpl->begin() + n);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::const_pointer
#else
CharT const*
#endif
basic_slim_string< CharT, TraitsT >::data() const { return m_pImpl->begin(); }

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::size() const { return m_pImpl->size(); }

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::const_iterator
#else
CharT const*
#endif
basic_slim_string< CharT, TraitsT >::begin() const { return m_pImpl->begin(); }
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::const_iterator
#else
CharT const*
#endif
basic_slim_string< CharT, TraitsT >::end() const { return m_pImpl->end(); }

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::copy(pointer s, size_type n, size_type pos) const
{
    const size_type size = m_pImpl->size();
    if (pos <= size)
    {
        const size_type len = (std::min)(n, size - pos);
        traits_type::copy(s, m_pImpl->begin() + pos, len);
        return len;
    }
    else
    {
        boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::copy: the position is out of range"));
    }
}

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find(basic_slim_string const& that, size_type pos) const
{
    return m_pImpl->find(that.m_pImpl->begin(), that.m_pImpl->size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find(string_type const& s, size_type pos) const
{
    return m_pImpl->find(s.data(), s.size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find(const_pointer s, size_type pos) const
{
    return m_pImpl->find(s, traits_type::length(s), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find(const_pointer s, size_type pos, size_type n) const
{
    return m_pImpl->find(s, n, pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find(char_type c, size_type pos) const
{
    return m_pImpl->find(c, pos);
}

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::rfind(basic_slim_string const& that, size_type pos) const
{
    return m_pImpl->rfind(that.m_pImpl->begin(), that.m_pImpl->size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::rfind(string_type const& s, size_type pos) const
{
    return m_pImpl->rfind(s.data(), s.size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::rfind(const_pointer s, size_type pos) const
{
    return m_pImpl->rfind(s, traits_type::length(s), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::rfind(const_pointer s, size_type pos, size_type n) const
{
    return m_pImpl->rfind(s, n, pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::rfind(char_type c, size_type pos) const
{
    return m_pImpl->rfind(c, pos);
}

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_of(basic_slim_string const& that, size_type pos) const
{
    return m_pImpl->find_first_of(that.m_pImpl->begin(), that.m_pImpl->size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_of(string_type const& s, size_type pos) const
{
    return m_pImpl->find_first_of(s.data(), s.size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_of(const_pointer s, size_type pos) const
{
    return m_pImpl->find_first_of(s, traits_type::length(s), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_of(const_pointer s, size_type pos, size_type n) const
{
    return m_pImpl->find_first_of(s, n, pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_of(char_type c, size_type pos) const
{
    return m_pImpl->find(c, pos);
}

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_of(basic_slim_string const& that, size_type pos) const
{
    return m_pImpl->find_last_of(that.m_pImpl->begin(), that.m_pImpl->size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_of(string_type const& s, size_type pos) const
{
    return m_pImpl->find_last_of(s.data(), s.size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_of(const_pointer s, size_type pos) const
{
    return m_pImpl->find_last_of(s, traits_type::length(s), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_of(const_pointer s, size_type pos, size_type n) const
{
    return m_pImpl->find_last_of(s, n, pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_of(char_type c, size_type pos) const
{
    return m_pImpl->rfind(c, pos);
}

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_not_of(basic_slim_string const& that, size_type pos) const
{
    return m_pImpl->find_first_not_of(that.m_pImpl->begin(), that.m_pImpl->size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_not_of(string_type const& s, size_type pos) const
{
    return m_pImpl->find_first_not_of(s.data(), s.size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_not_of(const_pointer s, size_type pos) const
{
    return m_pImpl->find_first_not_of(s, traits_type::length(s), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_not_of(const_pointer s, size_type pos, size_type n) const
{
    return m_pImpl->find_first_not_of(s, n, pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_first_not_of(char_type c, size_type pos) const
{
    return m_pImpl->find_first_not_of(c, pos);
}

template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_not_of(basic_slim_string const& that, size_type pos) const
{
    return m_pImpl->find_last_not_of(that.m_pImpl->begin(), that.m_pImpl->size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_not_of(string_type const& s, size_type pos) const
{
    return m_pImpl->find_last_not_of(s.data(), s.size(), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_not_of(const_pointer s, size_type pos) const
{
    return m_pImpl->find_last_not_of(s, traits_type::length(s), pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_not_of(const_pointer s, size_type pos, size_type n) const
{
    return m_pImpl->find_last_not_of(s, n, pos);
}
template< typename CharT, typename TraitsT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_slim_string< CharT, TraitsT >::size_type
#else
std::size_t
#endif
basic_slim_string< CharT, TraitsT >::find_last_not_of(char_type c, size_type pos) const
{
    return m_pImpl->find_last_not_of(c, pos);
}

template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(basic_slim_string const& that) const
{
    if (m_pImpl == that.m_pImpl)
        return 0;
    else
        return traits_type::compare(m_pImpl->begin(), that.m_pImpl->begin(), m_pImpl->size() + 1);
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(string_type const& s) const
{
    return traits_type::compare(m_pImpl->begin(), s.c_str(), m_pImpl->size() + 1);
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(size_type pos1, size_type n1, basic_slim_string const& that) const
{
    return m_pImpl->compare(pos1, n1, that.m_pImpl->begin(), that.m_pImpl->size());
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(size_type pos1, size_type n1, string_type const& s) const
{
    return m_pImpl->compare(pos1, n1, s.data(), s.size());
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(
    size_type pos1, size_type n1, basic_slim_string const& that, size_type pos2, size_type n2) const
{
    const size_type that_size = that.m_pImpl->size();
    if (that_size >= pos2)
    {
        return m_pImpl->compare(pos1, n1, that.m_pImpl->begin() + pos2, (std::min)(n2, that_size - pos2));
    }
    else
    {
        boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::compare: the position is out of range"));
    }
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(size_type pos1, size_type n1, string_type const& s, size_type pos2, size_type n2) const
{
    const size_type that_size = s.size();
    if (that_size >= pos2)
    {
        return m_pImpl->compare(pos1, n1, s.data() + pos2, (std::min)(n2, that_size - pos2));
    }
    else
    {
        boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::compare: the position is out of range"));
    }
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(const_pointer s) const
{
    return traits_type::compare(m_pImpl->begin(), s, m_pImpl->size() + 1);
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(const_pointer s, size_type n2) const
{
    return m_pImpl->compare(0, m_pImpl->size(), s, n2);
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(size_type pos1, size_type n1, const_pointer s) const
{
    // A bit more optimal version to avoid traits_type::length cycle
    const size_type size = m_pImpl->size();
    if (size >= pos1)
    {
        register const int res = traits_type::compare(m_pImpl->begin() + pos1, s, (std::min)(n1, size - pos1));
        if (res == 0)
            return (s[n1] == 0 ? 0 : -1);
        else
            return res;
    }
    else
    {
        boost::log::aux::throw_exception(std::out_of_range("basic_slim_string::compare: the position is out of range"));
    }
}
template< typename CharT, typename TraitsT >
int basic_slim_string< CharT, TraitsT >::compare(size_type pos1, size_type n1, const_pointer s, size_type n2) const
{
    return m_pImpl->compare(pos1, n1, s, n2);
}

#ifdef BOOST_LOG_USE_CHAR
template class basic_slim_string< char, std::char_traits< char > >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class basic_slim_string< wchar_t, std::char_traits< wchar_t > >;
#endif

} // namespace log

} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER
