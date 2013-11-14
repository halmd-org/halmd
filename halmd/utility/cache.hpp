/*
 * Copyright © 2013 Felix Höfling
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_UTILITY_CACHE_HPP
#define HALMD_UTILITY_CACHE_HPP

#include <cstdint>
#include <memory>
#include <type_traits>
#include <utility>

namespace halmd {

// forward declaration
template <typename T>
class cache_proxy;

/**
 * Cached value.
 *
 * This class provides a caching wrapper around an arbitrary value.
 * The state of the cache is maintained with an internal counter.
 * A write access to the value increments the counter, while a
 * read access does not change the counter.
 */
template <typename T = void>
class cache
{
private:
    /** roughly 10^39 cache writes before overflow */
    typedef std::uint64_t count_type;

public:
    /**
     * Construct cache.
     *
     * Forwards arguments to the value constructor without copying.
     */
    template <typename... Args>
    cache(Args&&... args)
      : value_(std::forward<Args>(args)...), count_(std::make_shared<count_type>(0)) {}

    /** deleted implicit copy constructor */
    cache(cache const&) = delete;
    /** deleted implicit assignment operator */
    cache& operator=(cache const&) = delete;
    /** define default move constructor */
    cache(cache&& c) = default;
    /** define default move assignment operator */
    cache& operator=(cache&& c) = default;

    /**
     * Returns const pointer to cached value for read access.
     */
    operator T const*() const
    {
        return &value_;
    }

    /**
     * Returns const pointer to cached value for read access.
     */
    T const* operator->() const
    {
        return &value_;
    }

    /**
     * Returns const reference to cached value for read access.
     */
    T const& operator*() const
    {
        return value_;
    }

private:
    friend class cache<void>;
    friend class cache_proxy<T>;
    friend class cache_proxy<T const>;

    /** cached value */
    T value_;
    /** cache write count */
    std::shared_ptr<count_type> count_;
};

/**
 * Observe cached value.
 *
 * This class may be used to check the state of a cached value.
 *
 * The observer supports transitive dependencies between cached values,
 * e.g. a chain of caches A → B → …. The equality operator of the cache
 * observer (of A) enforces the retrieval of the observed cache (A), by
 * calling a method that computes (from B) the value (A). This method may
 * in turn compare its own cache observer (of B) to the observed value (B),
 * which in turn enforces the retrieval of the observed cache (B), etc.
 */
template <>
class cache<void>
{
private:
    typedef std::uint64_t count_type;

public:
    /**
     * Default constructor.
     */
    cache() : initial_(-1UL) {}

    /**
     * Observe given cached value.
     */
    template <typename T>
    cache(cache<T> const& c) : count_(c.count_), initial_(*c.count_) {}

    /**
     * Returns true if cache is valid, false otherwise.
     */
    template <typename T>
    typename std::enable_if<!std::is_same<T, void>::value, bool>::type
    operator==(cache<T> const& c) const
    {
        return count_.lock() == c.count_ && initial_ == *c.count_;
    }

private:
    /** cache identity */
    std::weak_ptr<void> count_;
    /** cache write count at construction */
    count_type initial_;
};

/**
 * Returns true if cache is valid, false otherwise.
 */
template <typename T>
inline typename std::enable_if<!std::is_same<T, void>::value, bool>::type
operator==(cache<T> const& c, cache<> const& observer)
{
    return observer == c;
}

/**
 * Returns false if cache is valid, true otherwise.
 */
template <typename T>
inline typename std::enable_if<!std::is_same<T, void>::value, bool>::type
operator!=(cache<> const& observer, cache<T> const& c)
{
    return !(observer == c);
}

/**
 * Returns false if cache is valid, true otherwise.
 */
template <typename T>
inline typename std::enable_if<!std::is_same<T, void>::value, bool>::type
operator!=(cache<T> const& c, cache<> const& observer)
{
    return !(observer == c);
}

/**
 * Obtain write access to cached value.
 *
 * This class invalidates the cache.
 */
template <typename T>
class cache_proxy
{
public:
    typedef T& reference;
    typedef T* pointer;

    /**
     * Obtain write access to cached value.
     */
    cache_proxy(cache<T>& c) : value_(c.value_)
    {
        ++*c.count_;
    }

    /**
     * Returns reference to cached value.
     */
    reference operator*() const
    {
        return value_;
    }

    /**
     * Returns pointer to cached value.
     */
    pointer operator->() const
    {
        return &value_;
    }

    /** deleted implicit assignment operator */
    cache_proxy& operator=(cache_proxy const&) = delete;

private:
    /** reference to cached value */
    reference value_;
};

/**
 * Obtain read access to cached value.
 *
 * This class does not invalidate the cache.
 */
template <typename T>
class cache_proxy<T const>
{
public:
    typedef T const& reference;
    typedef T const* pointer;

    /**
     * Obtain read access to cached value.
     */
    cache_proxy(cache<T> const& c) : value_(c.value_) {}

    /**
     * Returns const reference to cached value.
     */
    reference operator*() const
    {
        return value_;
    }

    /**
     * Returns const pointer to cached value.
     */
    pointer operator->() const
    {
        return &value_;
    }

    /** deleted implicit assignment operator */
    cache_proxy& operator=(cache_proxy const&) = delete;

private:
    /** const reference to cached value */
    reference value_;
};

/**
 * Convenience function for write access to a cached value.
 *
 * Returns cache_proxy<T> by move constructor/assignment.
 *
 * Example:
 *
 *     cache<double> c(1.0);
 *     auto c_ = make_cache_mutable(c);
 *     *c_ = 2.0;
 *     assert(read_cache(c) == 2.0);
 */
template <typename T>
cache_proxy<T> make_cache_mutable(cache<T>& c) {
    return cache_proxy<T>(c);
}

/**
 * Convenience function for read access to a cached value.
 *
 * Returns const reference to cached value.
 *
 * Example:
 *
 *     const cache<double> c(1.0);
 *     auto const& c_ = read_cache(c);
 *     assert(c_ == 1.0);
 */
template <typename T>
T const& read_cache(cache<T> const& c) {
    return *cache_proxy<T const>(c);
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_CACHE_HPP */
