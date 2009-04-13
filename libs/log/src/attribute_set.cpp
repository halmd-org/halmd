/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   attribute_set.cpp
 * \author Andrey Semashev
 * \date   19.04.2007
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include <memory>
#include <functional>
#include <boost/intrusive/set.hpp>
#include <boost/intrusive/set_hook.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/list_hook.hpp>
#include <boost/intrusive/link_mode.hpp>
#include <boost/intrusive/derivation_value_traits.hpp>
#include <boost/utility/addressof.hpp>
#include <boost/log/attributes/attribute_set.hpp>
#include "light_key.hpp"

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace {

    //! A list-like container with ability to search for the value by the key
    template< typename ValueT, typename NodeT >
    class ordered_list :
        private std::allocator< ValueT >
    {
    public:
        typedef NodeT node_type;
        typedef std::allocator< ValueT > allocator_type;
        typedef typename allocator_type::value_type value_type;
        typedef typename allocator_type::reference reference;
        typedef typename allocator_type::const_reference const_reference;
        typedef typename allocator_type::pointer pointer;
        typedef typename allocator_type::const_pointer const_pointer;
        typedef typename allocator_type::size_type size_type;
        typedef typename allocator_type::difference_type difference_type;

    private:
        //! A simple functor that destroys the node and deallocates the memory
        struct disposer :
            public std::unary_function< pointer, void >
        {
            explicit disposer(allocator_type* pAlloc) : m_pAlloc(pAlloc) {}
            void operator() (pointer p) const
            {
                m_pAlloc->destroy(p);
                m_pAlloc->deallocate(p, 1);
            }

        private:
            allocator_type* m_pAlloc;
        };

        //! An ordering functor to allow the comparison of the nodes and keys
        struct templated_less
        {
            typedef bool result_type;

            template< typename T >
            bool operator() (const_reference left, T const& right) const
            {
                return (left.compare(right) < 0);
            }
            template< typename T >
            bool operator() (T const& left, const_reference right) const
            {
                return (right.compare(left) > 0);
            }
        };

        //! Node base class traits for the intrusive list
        struct node_traits
        {
            typedef node_type node;
            typedef node* node_ptr;
            typedef node const* const_node_ptr;
            static node* get_next(const node* n) { return n->m_pNext; }  
            static void set_next(node* n, node* next) { n->m_pNext = next; }  
            static node* get_previous(const node* n) { return n->m_pPrev; }  
            static void set_previous(node* n, node* prev) { n->m_pPrev = prev; }  
        };

        //! Contained node traits for the intrusive list
        typedef intrusive::derivation_value_traits< value_type, node_traits, intrusive::safe_link > value_traits;

        //! The intrusive list of nodes to fasten iteration
        typedef intrusive::list<
            value_type,
            intrusive::value_traits< value_traits >,
            intrusive::constant_time_size< false >
        > sequence_container;

        //! Ordered unique index to implement lookup
        typedef intrusive::set<
            value_type,
            intrusive::constant_time_size< true >
        > index_container;

    public:
        typedef typename sequence_container::iterator iterator;
        typedef typename sequence_container::const_iterator const_iterator;

    private:
        //! List of nodes. It is kept in the ordered state.
        sequence_container m_Sequence;
        //! Lookup index
        index_container m_Index;

    public:
        //! Default constructor
        ordered_list() {}
        //! Copy constructor
        ordered_list(ordered_list const& that) : allocator_type(static_cast< allocator_type const& >(that))
        {
            try
            {
                insert(that.begin(), that.end());
            }
            catch (...)
            {
                // If something happens we have to explicitly deallocate
                clear();
                throw;
            }
        }
        //! Destructor
        ~ordered_list() { clear(); }

        //! Assignment
        ordered_list& operator= (ordered_list const& that)
        {
            if (this != &that)
            {
                ordered_list tmp(that);
                swap(tmp);
            }
            return *this;
        }

        //! The method checks if the container is empty
        bool empty() const { return m_Sequence.empty(); }
        //! Returns the number of elements in the container
        size_type size() const { return m_Index.size(); }
        //! Swaps two instances to the container
        void swap(ordered_list& that)
        {
            m_Sequence.swap(that.m_Sequence);
            m_Index.swap(that.m_Index);
        }

        //  Iterator acquirement
        iterator begin() { return m_Sequence.begin(); }
        iterator end() { return m_Sequence.end(); }
        const_iterator begin() const { return m_Sequence.begin(); }
        const_iterator end() const { return m_Sequence.end(); }

        //! Clears the container
        void clear()
        {
            m_Index.clear();
            m_Sequence.clear_and_dispose(disposer(this));
        }

        //! Removes the element from the container
        void erase(iterator it)
        {
            m_Index.erase(m_Index.iterator_to(*it));
            m_Sequence.erase_and_dispose(it, disposer(this));
        }
        //! Removes the element from the container
        void erase(pointer p)
        {
            erase(m_Sequence.iterator_to(*p));
        }
        //! Removes the range of elements from the container
        void erase(iterator b, iterator e)
        {
            while (b != e)
                erase(b++);
        }
        //! Removes the element with the specified key from the container
        template< typename KeyT >
        size_type erase(KeyT const& key)
        {
            iterator it = find(key);
            if (it != end())
            {
                erase(it);
                return size_type(1);
            }
            else
                return 0;
        }

        //! Inserts the element into the container, if there is no other element with equivalent key
        std::pair< iterator, bool > insert(const_reference val)
        {
            typedef typename index_container::iterator index_iterator;
            index_iterator it = m_Index.lower_bound(val);
            pointer p = NULL;
            if (it == m_Index.end())
            {
                // No such element is in the container, ok to insert in the end
                p = construct_element(val);
                m_Sequence.push_back(*p);
            }
            else if (!it->is_equivalent(val))
            {
                // No such element is in the container, ok to insert in the middle
                p = construct_element(val);
                m_Sequence.insert(m_Sequence.iterator_to(*it), *p);
            }
            else
                return std::make_pair(m_Sequence.iterator_to(*it), false);

            it = m_Index.insert(it, *p);
            return std::make_pair(m_Sequence.iterator_to(*it), true);
        }
        //! Inserts the range of values into the container
        template< typename IteratorT >
        void insert(IteratorT b, IteratorT e)
        {
            for (; b != e; ++b)
                insert(*b);
        }

        //! Searches for the element with an equivalent key
        template< typename KeyT >
        iterator find(KeyT const& key)
        {
            typedef typename index_container::iterator index_iterator;
            index_iterator it = m_Index.find(key, templated_less());
            if (it != m_Index.end())
                return m_Sequence.iterator_to(*it);
            else
                return end();
        }
        //! Searches for the element with an equivalent key
        template< typename KeyT >
        const_iterator find(KeyT const& key) const
        {
            return const_iterator(const_cast< ordered_list* >(this)->find(key));
        }

    private:
        //! Allocates memory and constructs a copy of the element
        pointer construct_element(const_reference val)
        {
            pointer p = allocator_type::allocate(1);
            try
            {
                allocator_type::construct(p, val);
            }
            catch (...)
            {
                allocator_type::deallocate(p, 1);
                throw;
            }

            // These two won't throw
            sequence_container::node_algorithms::init(p);
            index_container::node_algorithms::init(p);

            return p;
        }
    };

} // namespace

//! Attribute set implementation
template< typename CharT >
struct basic_attribute_set< CharT >::implementation
{
public:
    //! A light key compound type
    typedef aux::light_key< char_type, size_type > light_key_type;
    
    //! The container node
    struct node :
        public node_base,
        public intrusive::set_base_hook< intrusive::link_mode< intrusive::safe_link > >
    {
        node(key_type const& key, mapped_type const& data) : node_base(key, data) {}

        bool operator< (node const& that) const
        {
            return (this->m_Value.first.compare(that.m_Value.first.c_str(), that.m_Value.first.size()) < 0);
        }
        int compare(light_key_type const& that) const
        {
            return this->m_Value.first.compare(that.pKey, that.KeyLen);
        }
        int compare(key_type const& that) const
        {
            return this->m_Value.first.compare(that.c_str(), that.size());
        }
        bool is_equivalent(node const& that) const
        {
            return (this->m_Value.first == that.m_Value.first);
        }
    };

    //! Node container type
    typedef ordered_list< node, node_base > node_container;

public:
    //! Node container
    node_container Nodes;
};

//! Default constructor
template< typename CharT >
basic_attribute_set< CharT >::basic_attribute_set() : m_pImpl(new implementation())
{
}

//! Copy constructor
template< typename CharT >
basic_attribute_set< CharT >::basic_attribute_set(basic_attribute_set const& that)
    : m_pImpl(new implementation(*that.m_pImpl))
{
}

//! Destructor
template< typename CharT >
basic_attribute_set< CharT >::~basic_attribute_set()
{
    delete m_pImpl;
}

//! Assignment
template< typename CharT >
basic_attribute_set< CharT >& basic_attribute_set< CharT >::operator= (basic_attribute_set const& that)
{
    if (this != &that)
    {
        basic_attribute_set tmp(that);
        swap(tmp);
    }
    return *this;
}

//  Iterator generators
template< typename CharT >
typename basic_attribute_set< CharT >::iterator basic_attribute_set< CharT >::begin()
{
    return iterator(m_pImpl->Nodes.begin().pointed_node());
}
template< typename CharT >
typename basic_attribute_set< CharT >::iterator basic_attribute_set< CharT >::end()
{
    return iterator(m_pImpl->Nodes.end().pointed_node());
}
template< typename CharT >
typename basic_attribute_set< CharT >::const_iterator basic_attribute_set< CharT >::begin() const
{
    return const_iterator(m_pImpl->Nodes.begin().pointed_node());
}
template< typename CharT >
typename basic_attribute_set< CharT >::const_iterator basic_attribute_set< CharT >::end() const
{
    return const_iterator(m_pImpl->Nodes.end().pointed_node());
}

//! The method returns number of elements in the container
template< typename CharT >
typename basic_attribute_set< CharT >::size_type basic_attribute_set< CharT >::size() const
{
    return m_pImpl->Nodes.size();
}

//! Insertion method
template< typename CharT >
std::pair< typename basic_attribute_set< CharT >::iterator, bool >
basic_attribute_set< CharT >::insert(key_type const& key, mapped_type const& data)
{
    typename implementation::node n(key, data);
    std::pair<
        typename implementation::node_container::iterator,
        bool
    > insertion_result = m_pImpl->Nodes.insert(n);

    return std::make_pair(iterator(insertion_result.first.pointed_node()), insertion_result.second);
}

//! The method erases all attributes with the specified name
template< typename CharT >
typename basic_attribute_set< CharT >::size_type
basic_attribute_set< CharT >::erase(key_type const& key)
{
    return m_pImpl->Nodes.erase(key);
}

//! The method erases the specified attribute
template< typename CharT >
void basic_attribute_set< CharT >::erase(iterator it)
{
    m_pImpl->Nodes.erase(static_cast< typename implementation::node* >(it.m_pNode));
}
//! The method erases all attributes within the specified range
template< typename CharT >
void basic_attribute_set< CharT >::erase(iterator begin, iterator end)
{
    while (begin != end)
        m_pImpl->Nodes.erase(static_cast< typename implementation::node* >((begin++).m_pNode));
}

//! The method clears the container
template< typename CharT >
void basic_attribute_set< CharT >::clear()
{
    m_pImpl->Nodes.clear();
}

//! Internal lookup implementation
template< typename CharT >
typename basic_attribute_set< CharT >::iterator
basic_attribute_set< CharT >::find_impl(const char_type* key, size_type len)
{
    return iterator(m_pImpl->Nodes.find(typename implementation::light_key_type(key, len)).pointed_node());
}

#ifdef BOOST_LOG_USE_CHAR
template class basic_attribute_set< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class basic_attribute_set< wchar_t >;
#endif

} // namespace log

} // namespace boost
