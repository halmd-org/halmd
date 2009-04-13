/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   attribute_values_view.cpp
 * \author Andrey Semashev
 * \date   19.04.2007
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <new>
#include <boost/assert.hpp>
#include <boost/log/attributes/attribute_values_view.hpp>
#ifndef BOOST_LOG_NO_THREADS
#include <boost/detail/atomic_count.hpp>
#endif // BOOST_LOG_NO_THREADS
#include "light_key.hpp"
#include "alignment_gap_between.hpp"

namespace boost {

namespace BOOST_LOG_NAMESPACE {

//! Container implementation
template< typename CharT >
struct basic_attribute_values_view< CharT >::implementation
{
public:
    //! Light generic key type
    typedef aux::light_key< char_type, size_type > light_key_type;

private:
    //! Ordering predicate to support light key in lookup
    struct light_key_less
    {
        typedef bool result_type;
        bool operator() (light_key_type const& left, node const& right) const
        {
            return (right.m_Value.first.compare(left.pKey, left.KeyLen) > 0);
        }
        bool operator() (node const& left, light_key_type const& right) const
        {
            return (left.m_Value.first.compare(right.pKey, right.KeyLen) < 0);
        }
    };

private:
    //! Reference counter
#ifndef BOOST_LOG_NO_THREADS
    detail::atomic_count m_RefCount;
#else
    unsigned long m_RefCount;
#endif // BOOST_LOG_NO_THREADS

    //! Points to the first element in the container
    node* m_pBegin;
    //! Points after the last element of the container
    node* m_pEnd;
    //! Points after the end of storage
    node* m_pEOS;

public:
    //! Constructor
    implementation(node* b, size_type n)
        : m_RefCount(1), m_pBegin(b), m_pEnd(b), m_pEOS(b + n)
    {
    }
    //! Destructor
    ~implementation()
    {
        for (node* p = m_pBegin, *e = m_pEnd; p != e; ++p)
            p->~node();
    }

    //! Returns the pointer to the first element
    node* begin() { return m_pBegin; }
    //! Returns the pointer after the last element
    node* end() { return m_pEnd; }

    //! Returns the number of elements in the container
    size_type size() const { return size_type(m_pEnd - m_pBegin); }
    //! Returns true if there are no elements in the container
    size_type capacity() const { return size_type(m_pEOS - m_pBegin); }

    //! Constructs an element at the end of the container. Does not check the storage or reallocate.
    void push_back(key_type const& key, attribute* attr)
    {
        BOOST_ASSERT(m_pEnd != m_pEOS);
        new (m_pEnd) node(key, attr);
        ++m_pEnd;
    }

    //! Looks for the element with an equivalent key
    node* find(light_key_type const& key) const
    {
        register node* b = m_pBegin;
        register node* e = m_pEnd;
        while (b != e)
        {
            register node* p = b + ((e - b) >> 1);
            register const int cmp = p->m_Value.first.compare(key.pKey, key.KeyLen);
            if (cmp == 0)
                return p;
            else if (cmp < 0)
                b = p + 1;
            else
                e = p;
        }

        return m_pEnd;
    }

    //! Increments the reference counter
    void add_ref() { ++m_RefCount; }
    //! Decrements the reference counter
    unsigned long release() { return static_cast< unsigned long >(--m_RefCount); }

    //! Constructs elements at the end of the view that correspond to the elements of the specified sequence
    void adopt_nodes(
        typename attribute_set_type::const_iterator& it,
        typename attribute_set_type::const_iterator end)
    {
        for (; it != end; ++it)
            push_back(it->first, it->second.get());
    }
    //! Constructs elements at the end of the view that correspond to the elements of the specified sequences.
    //! The function ensures the order of the created nodes and discards elements from [it2, end2) that have
    //! keys equivalent to the corresponding elements in [it1, end1).
    void adopt_nodes(
        typename attribute_set_type::const_iterator& it1,
        typename attribute_set_type::const_iterator end1,
        typename attribute_set_type::const_iterator& it2,
        typename attribute_set_type::const_iterator end2)
    {
        while (true)
        {
            unsigned int cond =
                (static_cast< unsigned int >(it1 != end1) << 1)
                | static_cast< unsigned int >(it2 != end2);

            switch (cond)
            {
            case 1:
                adopt_nodes(it2, end2);
                return;

            case 2:
                adopt_nodes(it1, end1);
                return;

            case 3:
            {
                register int cmp = it1->first.compare(it2->first);
                if (cmp > 0)
                {
                    push_back(it2->first, it2->second.get());
                    ++it2;
                }
                else
                {
                    push_back(it1->first, it1->second.get());
                    ++it1;

                    if (cmp == 0) // we skip equivalent-keyed elements in the lesser priority sets
                        ++it2;
                }
                break;
            }

            default:
                return;
            }
        }
    }

    //! Allocates the needed memory to accomodate n nodes and constructs implementation object with m_RefCount == 1.
    static implementation* allocate_and_add_ref(internal_allocator_type* pAlloc, size_type n)
    {
        enum { alignment_gap = aux::alignment_gap_between< implementation, node >::value };
        typename internal_allocator_type::pointer p =
            pAlloc->allocate(sizeof(implementation) + (n * sizeof(node)) + alignment_gap);
        new (p) implementation(reinterpret_cast< node* >(p + sizeof(implementation) + alignment_gap), n);
        return reinterpret_cast< implementation* >(p);
    }
    //! Releases the implementation object and if its ref. counter drops to 0, destroys it and deallocates memory.
    static void release_and_dispose(implementation* pImpl, internal_allocator_type* pAlloc)
    {
        if (pImpl->release() == 0)
        {
            enum { alignment_gap = aux::alignment_gap_between< implementation, node >::value };
            const size_type size =
                sizeof(implementation) + (pImpl->capacity() * sizeof(node)) + alignment_gap;
            pImpl->~implementation();
            pAlloc->deallocate(
                reinterpret_cast< typename internal_allocator_type::pointer >(pImpl),
                size);
        }
    }
};

//! The constructor adopts three attribute sets to the view
template< typename CharT >
basic_attribute_values_view< CharT >::basic_attribute_values_view(
    attribute_set_type const& source_attrs,
    attribute_set_type const& thread_attrs,
    attribute_set_type const& global_attrs)
{
    // Allocate the implementation container
    size_type max_size = source_attrs.size() + thread_attrs.size() + global_attrs.size();
    m_pImpl = implementation::allocate_and_add_ref(this, max_size);

    if (max_size > 0)
    {
        // Compose the view. All copies performed bellow don't throw, so we are safe now.
        typename attribute_set_type::const_iterator
            src = source_attrs.begin(),
            esrc = source_attrs.end(),
            trd = thread_attrs.begin(),
            etrd = thread_attrs.end(),
            glb = global_attrs.begin(),
            eglb = global_attrs.end();
    
        while (true)
        {
            unsigned int cond =
                (static_cast< unsigned int >(src != esrc) << 2)
                | (static_cast< unsigned int >(trd != etrd) << 1)
                | static_cast< unsigned int >(glb != eglb);

            switch (cond)
            {
            case 1:
                // Only global attributes left
                m_pImpl->adopt_nodes(glb, eglb);
                return;

            case 2:
                // Only thread attributes left
                m_pImpl->adopt_nodes(trd, etrd);
                return;

            case 3:
                // Only thread and global attributes left
                m_pImpl->adopt_nodes(trd, etrd, glb, eglb);
                return;

            case 4:
                // Only source attributes left
                m_pImpl->adopt_nodes(src, esrc);
                return;

            case 5:
                // Only source and global attributes left
                m_pImpl->adopt_nodes(src, esrc, glb, eglb);
                return;

            case 6:
                // Only source and thread attributes left
                m_pImpl->adopt_nodes(src, esrc, trd, etrd);
                return;

            case 7:
            {
                // Source, thread and global attributes left
                // Select the least-keyed attribute
                register typename attribute_set_type::const_iterator* pIt = &src;
                register int cmp = (*pIt)->first.compare(trd->first);
                if (cmp > 0)
                    pIt = &trd;
                else if (cmp == 0) // we skip equivalent-keyed elements in the lesser priority sets
                    ++trd;
                cmp = (*pIt)->first.compare(glb->first);
                if (cmp > 0)
                    pIt = &glb;
                else if (cmp == 0) // we skip equivalent-keyed elements in the lesser priority sets
                    ++glb;

                // Adopt the attribute
                m_pImpl->push_back((*pIt)->first, (*pIt)->second.get());
                ++(*pIt);

                break;
            }

            default:
                // No attributes left to adopt
                return;
            }
        }
    }
}

//! Copy constructor
template< typename CharT >
basic_attribute_values_view< CharT >::basic_attribute_values_view(basic_attribute_values_view const& that)
    : internal_allocator_type(static_cast< internal_allocator_type const& >(that)), m_pImpl(that.m_pImpl)
{
    m_pImpl->add_ref();
}

//! Destructor
template< typename CharT >
basic_attribute_values_view< CharT >::~basic_attribute_values_view()
{
    implementation::release_and_dispose(m_pImpl, this);
}

//! Assignment
template< typename CharT >
basic_attribute_values_view< CharT >& basic_attribute_values_view< CharT >::operator= (basic_attribute_values_view const& that)
{
    if (this != &that)
    {
        basic_attribute_values_view tmp(that);
        swap(tmp);
    }
    return *this;
}

//  Iterator generators
template< typename CharT >
typename basic_attribute_values_view< CharT >::const_iterator
basic_attribute_values_view< CharT >::begin() const
{
    return const_iterator(m_pImpl->begin());
}

template< typename CharT >
typename basic_attribute_values_view< CharT >::const_iterator
basic_attribute_values_view< CharT >::end() const
{
    return const_iterator(m_pImpl->end());
}

//! The method returns number of elements in the container
template< typename CharT >
typename basic_attribute_values_view< CharT >::size_type basic_attribute_values_view< CharT >::size() const
{
    return m_pImpl->size();
}

//! Internal lookup implementation
template< typename CharT >
typename basic_attribute_values_view< CharT >::const_iterator
basic_attribute_values_view< CharT >::find_impl(const char_type* key, size_type len) const
{
    return const_iterator(m_pImpl->find(typename implementation::light_key_type(key, len)));
}

//! The method acquires values of all adopted attributes. Users don't need to call it, since will always get an already frozen view.
template< typename CharT >
void basic_attribute_values_view< CharT >::freeze()
{
    for (const_iterator it = begin(), e = end(); it != e; ++it)
        it.freeze_element();
}

//! Explicitly instantiate container implementation
#ifdef BOOST_LOG_USE_CHAR
template class basic_attribute_values_view< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class basic_attribute_values_view< wchar_t >;
#endif

} // namespace log

} // namespace boost
