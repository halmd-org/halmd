/*!
 * (C) 2008 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   new_shared.hpp
 * \author Andrey Semashev
 * \date   08.05.2008
 * 
 * \brief  This header implements a helper function new_shared that allocates an
 *         object and returns a shared_ptr to it. It may perform better than
 *         doing so with explicit call to operator new and shared_ptr constructor.
 *         The new_shared function may be taken as a replacement for the core operator new.
 */

#ifndef BOOST_LOG_NEW_SHARED_HPP_INCLUDED_
#define BOOST_LOG_NEW_SHARED_HPP_INCLUDED_

#include <cstddef>
#include <new>
#include <memory>
#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/aligned_storage.hpp>
#include <boost/type_traits/alignment_of.hpp>

#if !(defined(BOOST_HAS_VARIADIC_TMPL) && defined(BOOST_HAS_RVALUE_REFS))

#include <boost/preprocessor/comma_if.hpp>
#include <boost/preprocessor/arithmetic/dec.hpp>
#include <boost/preprocessor/comparison/greater.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

#ifndef BOOST_LOG_NEW_SHARED_MAX_CTOR_ARGS
#define BOOST_LOG_NEW_SHARED_MAX_CTOR_ARGS 10
#endif // BOOST_LOG_NEW_SHARED_MAX_CTOR_ARGS

#endif // !(defined(BOOST_HAS_VARIADIC_TMPL) && defined(BOOST_HAS_RVALUE_REFS))

#include <boost/log/detail/prologue.hpp>

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace aux {

//! Default deleter for the pointer returned from new_shared
struct new_shared_deleter
{
    typedef void result_type;
    template< typename T >
    void operator() (T* p) const
    {
        p->~T();
    }
};

//! Wrapper allocator for the shared_ptr counter objects
template< typename T, typename StorageT >
struct new_shared_allocator :
    private StorageT::allocator_type
{
    template< typename, typename >
    friend struct new_shared_allocator;

private:
    typedef StorageT storage_type;
    typedef typename storage_type::allocator_type base_type;

public:
    typedef T value_type;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef T& reference;
    typedef T const& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    template< typename U >
    struct rebind
    {
        // Check that the storage is capable to hold the new counter value
        BOOST_STATIC_ASSERT(sizeof(U) <= sizeof(typename storage_type::sp_counter_type));
        BOOST_STATIC_ASSERT(alignment_of< U >::value <= alignment_of< typename storage_type::sp_counter_type >::value);
        typedef new_shared_allocator< U, storage_type > other;
    };

public:
    new_shared_allocator(storage_type* pStorage, base_type const& a) :
        base_type(a),
        m_pStorage(pStorage)
    {
    }
    new_shared_allocator(new_shared_allocator const& that) :
        base_type(static_cast< base_type const& >(that)),
        m_pStorage(that.m_pStorage)
    {
        that.m_pStorage = 0;
    }
    template< typename U >
    new_shared_allocator(new_shared_allocator< U, storage_type > const& that) :
        base_type(static_cast< base_type const& >(that)),
        m_pStorage(that.m_pStorage)
    {
        that.m_pStorage = 0;
    }
    ~new_shared_allocator()
    {
        // Need to check m_pStorage since deallocating NULL of size 1 is illegal
        if (m_pStorage)
            base_type::deallocate(m_pStorage, 1);
    }

    pointer address(reference r) const { return &r; }
    const_pointer address(const_reference r) const { return &r; }
    size_type max_size() const { return base_type::max_size(); }

    pointer allocate(size_type n, const void* = 0)
    {
        if (m_pStorage && n == 1)
        {
            pointer p = reinterpret_cast< pointer >(&m_pStorage->m_CounterStorage);
            m_pStorage = 0;
            return p;
        }
        else
            throw std::bad_alloc();
    }
    void deallocate(pointer p, size_type n)
    {
        if (p)
        {
            base_type::deallocate(reinterpret_cast< typename base_type::pointer >(
                reinterpret_cast< char* >(p) - offsetof(storage_type, m_CounterStorage)), 1);
        }
    }

    void construct(pointer p, const_reference val) { new(p) value_type(val); }
    void destroy(pointer p) { p->~value_type(); }

private:
    mutable storage_type* m_pStorage;
};

//! The storage for both the shared_ptr counter object and the pointed object itself (it must be POD)
template<
    typename T,
    typename DeleterT,
    typename AllocatorT
>
struct new_shared_storage
{
    typedef typename AllocatorT::BOOST_NESTED_TEMPLATE rebind< new_shared_storage >::other allocator_type;
    typedef boost::detail::sp_counted_impl_pda< T*, DeleterT, new_shared_allocator< T, new_shared_storage > > sp_counter_type;
    typedef DeleterT deleter_type;

    //! The storage for the counter (it must be POD)
    typename boost::aligned_storage<
        sizeof(sp_counter_type),
        alignment_of< sp_counter_type >::value
    >::type m_CounterStorage;

    //! The storage for the pointed object (it must be POD)
    typename boost::aligned_storage<
        sizeof(T),
        alignment_of< T >::value
    >::type m_ValueStorage;
};


template< typename DeleterT >
struct use_deleter_param
{
    DeleterT value;
    explicit use_deleter_param(DeleterT const& d) : value(d) {}
};

template< typename AllocatorT >
struct use_allocator_param
{
    AllocatorT value;
    explicit use_allocator_param(AllocatorT const& a) : value(a) {}
};


template< typename DeleterT >
inline use_deleter_param< DeleterT > use_deleter(DeleterT const& d)
{
    return use_deleter_param< DeleterT >(d);
}

template< typename AllocatorT >
inline use_allocator_param< AllocatorT > use_allocator(AllocatorT const& a)
{
    return use_allocator_param< AllocatorT >(a);
}

#if (defined(BOOST_HAS_VARIADIC_TMPL) && defined(BOOST_HAS_RVALUE_REFS))

template<
    typename T,
    typename DeleterT,
    typename AllocatorT,
    typename... ArgsT
>
inline shared_ptr< T > new_shared_impl(DeleterT const& d, AllocatorT const& a, ArgsT&&... args)
{
    typedef new_shared_storage< T, DeleterT, AllocatorT > storage_type;
    typedef typename storage_type::allocator_type allocator_type;
    typedef typename storage_type::deleter_type deleter_type;

    allocator_type rebound_alloc(a);
    storage_type* pStorage = rebound_alloc.allocate(1);
    T* pObj = reinterpret_cast< T* >(&pStorage->m_ValueStorage);

    try
    {
        new(pObj) T(args...);
    }
    catch (...)
    {
        rebound_alloc.deallocate(pStorage, 1);
        throw;
    }

    // The new_shared_allocator won't throw
    // And just in case if shared_ptr constructor throws for some another reason
    // new_shared_allocator deallocates storage on its destructor.
    return shared_ptr< T >(
        pObj,
        d,
        new_shared_allocator< T, storage_type >(pStorage, rebound_alloc));
}

template< typename T >
inline shared_ptr< T > new_shared()
{
    new_shared_deleter d;
    std::allocator< void > a;
    return new_shared_impl< T >(d, a);
}

template< typename T, typename... ArgsT >
inline shared_ptr< T > new_shared(ArgsT&&... args)
{
    new_shared_deleter d;
    std::allocator< void > a;
    return new_shared_impl< T >(d, a, args...);
}

template< typename T, typename DeleterT, typename... ArgsT >
inline shared_ptr< T > new_shared(use_deleter_param< DeleterT >&& d, ArgsT&&... args)
{
    std::allocator< void > a;
    return new_shared_impl< T >(d.value, a, args...);
}

template< typename T, typename AllocatorT, typename... ArgsT >
inline shared_ptr< T > new_shared(use_allocator_param< AllocatorT >&& a, ArgsT&&... args)
{
    new_shared_deleter d;
    return new_shared_impl< T >(d, a.value, args...);
}

template< typename T, typename DeleterT, typename AllocatorT, typename... ArgsT >
inline shared_ptr< T > new_shared(
    use_deleter_param< DeleterT >&& d,
    use_allocator_param< AllocatorT >&& a,
    ArgsT&&... args)
{
    return new_shared_impl< T >(d.value, a.value, args...);
}

template< typename T, typename DeleterT, typename AllocatorT, typename... ArgsT >
inline shared_ptr< T > new_shared(
    use_allocator_param< AllocatorT >&& a,
    use_deleter_param< DeleterT >&& d,
    ArgsT&&... args)
{
    return new_shared_impl< T >(d.value, a.value, args...);
}

#else // (defined(BOOST_HAS_VARIADIC_TMPL) && defined(BOOST_HAS_RVALUE_REFS))

#define BOOST_PP_FILENAME_1 <boost/log/detail/new_shared.hpp>
#define BOOST_PP_ITERATION_LIMITS (0, BOOST_LOG_NEW_SHARED_MAX_CTOR_ARGS)
#include BOOST_PP_ITERATE()

#endif // (defined(BOOST_HAS_VARIADIC_TMPL) && defined(BOOST_HAS_RVALUE_REFS))

} // namespace aux

} // namespace log

} // namespace boost

#endif // BOOST_NEW_SHARED_HPP_INCLUDED_

#if !(defined(BOOST_HAS_VARIADIC_TMPL) && defined(BOOST_HAS_RVALUE_REFS)) && defined(BOOST_PP_IS_ITERATING)

template<
    typename T,
    typename DeleterT,
    typename AllocatorT BOOST_PP_COMMA_IF(BOOST_PP_ITERATION())
    BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), typename ArgT)
>
inline shared_ptr< T > new_shared_impl(
    DeleterT const& d,
    AllocatorT const& a BOOST_PP_COMMA_IF(BOOST_PP_ITERATION())
    BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_ITERATION(), ArgT, const& arg)
)
{
    typedef new_shared_storage< T, DeleterT, AllocatorT > storage_type;
    typedef typename storage_type::allocator_type allocator_type;
    typedef typename storage_type::deleter_type deleter_type;

    allocator_type rebound_alloc(a);
    storage_type* pStorage = rebound_alloc.allocate(1);
    T* pObj = reinterpret_cast< T* >(&pStorage->m_ValueStorage);

    try
    {
        new(pObj) T(BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), arg));
    }
    catch (...)
    {
        rebound_alloc.deallocate(pStorage, 1);
        throw;
    }

    // The new_shared_allocator won't throw
    // And just in case if shared_ptr constructor throws for some another reason
    // new_shared_allocator deallocates storage on its destructor.
    return shared_ptr< T >(
        pObj,
        d,
        new_shared_allocator< T, storage_type >(pStorage, rebound_alloc));
}

#if BOOST_PP_ITERATION() == 0

template< typename T >
inline shared_ptr< T > new_shared()
{
    new_shared_deleter d;
    std::allocator< void > a;
    return new_shared_impl< T >(d, a);
}

#elif BOOST_PP_ITERATION() == 1

template< typename T, typename ArgT0 >
inline shared_ptr< T > new_shared(ArgT0 const& arg0)
{
    new_shared_deleter d;
    std::allocator< void > a;
    return new_shared_impl< T >(d, a, arg0);
}

template< typename T, typename DeleterT >
inline shared_ptr< T > new_shared(use_deleter_param< DeleterT > const& d)
{
    std::allocator< void > a;
    return new_shared_impl< T >(d.value, a);
}

template< typename T, typename AllocatorT >
inline shared_ptr< T > new_shared(use_allocator_param< AllocatorT > const& a)
{
    new_shared_deleter d;
    return new_shared_impl< T >(d, a.value);
}

#else

template< typename T, BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), typename ArgT) >
inline shared_ptr< T > new_shared(BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_ITERATION(), ArgT, const& arg))
{
    new_shared_deleter d;
    std::allocator< void > a;
    return new_shared_impl< T >(d, a, BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), arg));
}

template< typename T, typename DeleterT, BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_ITERATION()), typename ArgT) >
inline shared_ptr< T > new_shared(
    use_deleter_param< DeleterT > const& d,
    BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_DEC(BOOST_PP_ITERATION()), ArgT, const& arg))
{
    std::allocator< void > a;
    return new_shared_impl< T >(d.value, a, BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_ITERATION()), arg));
}

template< typename T, typename AllocatorT, BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_ITERATION()), typename ArgT) >
inline shared_ptr< T > new_shared(
    use_allocator_param< AllocatorT > const& a,
    BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_DEC(BOOST_PP_ITERATION()), ArgT, const& arg))
{
    new_shared_deleter d;
    return new_shared_impl< T >(d, a.value, BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_ITERATION()), arg));
}

template<
    typename T,
    typename DeleterT,
    typename AllocatorT BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ITERATION(), 2))
    BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_ITERATION())), typename ArgT)
>
inline shared_ptr< T > new_shared(
    use_deleter_param< DeleterT > const& d,
    use_allocator_param< AllocatorT > const& a BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ITERATION(), 2))
    BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_ITERATION())), ArgT, const& arg))
{
    return new_shared_impl< T >(
        d.value,
        a.value BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ITERATION(), 2))
        BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_ITERATION())), arg));
}

template<
    typename T,
    typename DeleterT,
    typename AllocatorT BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ITERATION(), 2))
    BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_ITERATION())), typename ArgT)
>
inline shared_ptr< T > new_shared(
    use_allocator_param< AllocatorT > const& a,
    use_deleter_param< DeleterT > const& d BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ITERATION(), 2))
    BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_ITERATION())), ArgT, const& arg))
{
    return new_shared_impl< T >(
        d.value,
        a.value BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ITERATION(), 2))
        BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(BOOST_PP_DEC(BOOST_PP_ITERATION())), arg));
}

#endif // BOOST_PP_ITERATION() == n

#endif // !(defined(BOOST_HAS_VARIADIC_TMPL) && defined(BOOST_HAS_RVALUE_REFS)) && defined(BOOST_PP_IS_ITERATING)
