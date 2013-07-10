/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#ifndef HALMD_OBSERVABLES_DYNAMICS_CORRELATION_HPP
#define HALMD_OBSERVABLES_DYNAMICS_CORRELATION_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

/**
 * Store input samples (phase space, density modes, ...) in a
 * coarse-grained block structure and provide a signal
 * on_correlate_block which correlation functions connect to.
 */
class correlation_base
{
public:
    correlation_base() {}
    virtual ~correlation_base() {}
    /** Lua bindings */
    static void luaopen(lua_State* L);

    /** compute correlations at the given coarse-graining level */
    virtual void compute(unsigned int level) = 0;
};

namespace detail {

/**
 * determine result_rank of generic TCF functor or provide a default
 */
template <typename T>
constexpr auto get_rank_impl(int) -> decltype(T::result_rank)
{
    return T::result_rank;
}

template <typename T>
constexpr int get_rank_impl(long)
{
    return 0;
}

template <typename T>
constexpr int get_rank()
{
    return get_rank_impl<T>(0); // 0 is of type 'int' which takes precedence over 'long'
}

/**
 * call result_shape() of generic TCF functor or provide a default
 */
template <typename T>
auto get_shape_impl(T const& obj, int) -> decltype(obj.result_shape())
{
    return obj.result_shape();
}

template <typename T>
unsigned int const* get_shape_impl(T const& obj, long)
{
    return nullptr;
}

template <typename T>
unsigned int const* get_shape(T const& obj)
{
    return get_shape_impl(obj, 0); // 0 is of type 'int' which takes precedence over 'long'
}

} // namespace detail

template <typename tcf_type>
class correlation
  : public correlation_base
{
    enum { rank_ = 2 + detail::get_rank<tcf_type>() };

public:
    typedef typename tcf_type::sample_type sample_type;
    typedef typename tcf_type::result_type result_type;

    typedef observables::samples::blocking_scheme<sample_type> block_sample_type;
    typedef boost::multi_array<accumulator<result_type>, rank_> block_result_type;
    typedef boost::multi_array<result_type, rank_> block_mean_type;
    typedef boost::multi_array<typename accumulator<result_type>::size_type, rank_> block_count_type;
    typedef logger logger_type;

    static void luaopen(lua_State* L);

    correlation(
        std::shared_ptr<tcf_type> tcf
      , std::shared_ptr<block_sample_type> block_sample1
      , std::shared_ptr<block_sample_type> block_sample2
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>("dynamics.correlation")
    );

    virtual void compute(unsigned int level);

    block_result_type const& result() const
    {
        return result_;
    }

    static std::function<block_mean_type const& ()>
    get_mean(std::shared_ptr<correlation> self);

    static std::function<block_mean_type const& ()>
    get_error(std::shared_ptr<correlation> self);

    static std::function<block_count_type const& ()>
    get_count(std::shared_ptr<correlation> self);

private:
    typedef correlation_base _Base;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        utility::profiler::accumulator_type tcf;
    };

    /** block structures holding the input data */
    std::shared_ptr<block_sample_type> block_sample1_;
    std::shared_ptr<block_sample_type> block_sample2_;
    /** functor performing the specific computation */
    std::shared_ptr<tcf_type> tcf_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;

    /** block structures holding accumulated result values */
    block_result_type result_;
    /** mean values */
    block_mean_type mean_;
    /** standard error of mean */
    block_mean_type error_;
    /** accumulator count */
    block_count_type count_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

template <typename tcf_type>
correlation<tcf_type>::correlation(
    std::shared_ptr<tcf_type> tcf
  , std::shared_ptr<block_sample_type> block_sample1
  , std::shared_ptr<block_sample_type> block_sample2
  , std::shared_ptr<logger_type> logger
)
  // dependency injection
  : block_sample1_(block_sample1)
  , block_sample2_(block_sample2)
  , tcf_(tcf)
  , logger_(logger)
{
    if (block_sample1_->count() != block_sample2_->count()
     || block_sample1_->block_size() != block_sample2_->block_size())
    {
        LOG_ERROR("blocking schemes have incompatible shapes");
        throw std::logic_error("failed to setup time correlation function.");
    }

    // construct shape of total result array from shape of correlation function
    boost::array<typename block_result_type::size_type, rank_> extents;
    extents[0] = block_sample1_->count();
    extents[1] = block_sample1_->block_size();
    auto shape = detail::get_shape(*tcf_);
    for (unsigned int i = 2; i < rank_; ++i) {
        extents[i] = shape[i - 2];
    }
    // memory allocation
    result_.resize(extents);
    mean_.resize(extents);
    error_.resize(extents);
    count_.resize(extents);
}

template <typename tcf_type>
void correlation<tcf_type>::compute(unsigned int level)
{
    LOG_TRACE("compute correlations at level " << level);

    // iterate over block and correlate the first entry of block_sample1_ (at
    // time t1) with all entries of block_sample2_ (at t1 + n * Δt), accumulate
    // result for each lag time
    auto out = result_[level].begin();
    auto first = block_sample1_->index(level).begin();
    auto const& block2 = block_sample2_->index(level);
    for (auto second = block2.begin(); second != block2.end(); ++second) {
        scoped_timer_type timer(runtime_.tcf);
        // call TCF functor which correlates the two samples and stores the
        // result in the output accumulator
        (*tcf_)(**first, **second, *out++);
    }
}

template <typename tcf_type>
std::function<typename correlation<tcf_type>::block_mean_type const& ()>
correlation<tcf_type>::get_mean(std::shared_ptr<correlation<tcf_type>> self)
{
    return [=]() -> block_mean_type const& {
        auto in  = self->result_.origin();
        auto out = self->mean_.origin();
        for (unsigned int i = 0; i < self->mean_.num_elements(); ++i) {
            *out++ = mean(*in++);
        }
        return self->mean_;
    };
}

template <typename tcf_type>
std::function<typename correlation<tcf_type>::block_mean_type const& ()>
correlation<tcf_type>::get_error(std::shared_ptr<correlation<tcf_type>> self)
{
    return [=]() -> block_mean_type const& {
        auto in  = self->result_.origin();
        auto out = self->error_.origin();
        for (unsigned int i = 0; i < self->error_.num_elements(); ++i) {
            *out++ = error_of_mean(*in++);
        }
        return self->error_;
    };
}

template <typename tcf_type>
std::function<typename correlation<tcf_type>::block_count_type const& ()>
correlation<tcf_type>::get_count(std::shared_ptr<correlation<tcf_type>> self)
{
    return [=]() -> block_count_type const& {
        auto in  = self->result_.origin();
        auto out = self->count_.origin();
        for (unsigned int i = 0; i < self->count_.num_elements(); ++i) {
            *out++ = count(*in++);
        }
        return self->count_;
    };
}

template <typename tcf_type>
void correlation<tcf_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<correlation, _Base>()
                    .property("mean", &correlation::get_mean)
                    .property("error", &correlation::get_error)
                    .property("count", &correlation::get_count)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("tcf", &runtime::tcf)
                    ]
                    .def_readonly("runtime", &correlation::runtime_)

              , def("correlation", &std::make_shared<correlation
                  , std::shared_ptr<tcf_type>
                  , std::shared_ptr<block_sample_type>
                  , std::shared_ptr<block_sample_type>
                  , std::shared_ptr<logger_type>
                >)
            ]
        ]
    ];
}

} // namespace dynamics
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_CORRELATION_HPP */
