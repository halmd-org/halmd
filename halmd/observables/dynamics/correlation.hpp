/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>

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

template <typename tcf_type>
class correlation
  : public correlation_base
{
public:
    typedef typename tcf_type::sample_type sample_type;
    typedef typename tcf_type::result_type result_type;
    typedef accumulator<result_type> accumulator_type;
    typedef observables::samples::blocking_scheme<sample_type> block_sample_type;
    typedef boost::multi_array<accumulator_type, 2> block_result_type;
    typedef boost::multi_array<typename accumulator_type::value_type, 2> block_mean_type;
    typedef boost::multi_array<typename accumulator_type::value_type, 2> block_error_type;
    typedef boost::multi_array<typename accumulator_type::size_type, 2> block_count_type;

    static void luaopen(lua_State* L, char const* scope);

    correlation(
        boost::shared_ptr<tcf_type> tcf
      , boost::shared_ptr<block_sample_type> block_sample
    );
    virtual ~correlation() {}

    virtual void compute(unsigned int level);
    block_mean_type const& mean();
    block_error_type const& error();
    block_count_type const& count();

private:
    typedef correlation_base _Base;

    /** block structure holding the input data */
    boost::shared_ptr<block_sample_type> block_sample_;
    /** functor performing the specific computation */
    boost::shared_ptr<tcf_type> tcf_;
    /** block structures holding accumulated result values */
    block_result_type result_;
    /** mean values */
    block_mean_type mean_;
    /** standard error of mean */
    block_error_type error_;
    /** accumulator count */
    block_count_type count_;
};

template <typename tcf_type>
correlation<tcf_type>::correlation(
    boost::shared_ptr<tcf_type> tcf
  , boost::shared_ptr<block_sample_type> block_sample
)
  // dependency injection
  : block_sample_(block_sample)
  , tcf_(tcf)
  // memory allocation
  , result_(boost::extents[block_sample->count()][block_sample->block_size()])
  , mean_(boost::extents[block_sample->count()][block_sample->block_size()])
  , error_(boost::extents[block_sample->count()][block_sample->block_size()])
  , count_(boost::extents[block_sample->count()][block_sample->block_size()])
{
}

template <typename tcf_type>
void correlation<tcf_type>::compute(unsigned int level)
{
    LOG_TRACE("[" << tcf_type::module_name() << "]: compute correlations at level " << level);

    typedef typename block_sample_type::block_type block_type;
    typedef typename block_type::const_iterator input_iterator;
    typedef typename block_result_type::reference::iterator output_iterator;

    // iterate over block and correlate the first entry (at time t1)
    // with all entries (at t1 + n * Δt), accumulate result for each lag time
    block_type const& block = block_sample_->index(level);
    input_iterator first = block.begin();
    output_iterator out = result_[level].begin();
    for (input_iterator second = first; second != block.end(); ++second) {
        // call TCF-specific compute routine and
        // store result in output accumulator
        (*out++)(tcf_->compute(*first, *second));
    }
}

template <typename tcf_type>
typename correlation<tcf_type>::block_mean_type const&
correlation<tcf_type>::mean()
{
    for (std::size_t i = 0; i < result_.shape()[0]; ++i) {
        for (std::size_t j = 0; j < result_.shape()[1]; ++j) {
            mean_[i][j] = detail::numeric::mean(result_[i][j]);
        }
    }
    return mean_;
}

template <typename tcf_type>
typename correlation<tcf_type>::block_error_type const&
correlation<tcf_type>::error()
{
    for (std::size_t i = 0; i < result_.shape()[0]; ++i) {
        for (std::size_t j = 0; j < result_.shape()[1]; ++j) {
            error_[i][j] = error_of_mean(result_[i][j]);
        }
    }
    return error_;
}

template <typename tcf_type>
typename correlation<tcf_type>::block_count_type const&
correlation<tcf_type>::count()
{
    for (std::size_t i = 0; i < result_.shape()[0]; ++i) {
        for (std::size_t j = 0; j < result_.shape()[1]; ++j) {
            count_[i][j] = detail::numeric::count(result_[i][j]);
        }
    }
    return count_;
}

template <typename tcf_type>
static char const* class_name_wrapper(correlation<tcf_type> const&)
{
    return tcf_type::class_name();
}

template <typename tcf_type>
static char const* module_name_wrapper(correlation<tcf_type> const&)
{
    return tcf_type::module_name();
}

template <typename tcf_type>
static char const* sample_name_wrapper(correlation<tcf_type> const&)
{
    return tcf_type::sample_type::class_name();
}

template <typename tcf_type>
void correlation<tcf_type>::luaopen(lua_State* L, char const* scope)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_(scope)
            [
                namespace_("dynamics")
                [
                    namespace_("correlation")
                    [
                        class_<correlation, boost::shared_ptr<_Base>, _Base>(tcf_type::class_name())
                            .def(constructor<
                                boost::shared_ptr<tcf_type>
                              , boost::shared_ptr<block_sample_type>
                            >())
                            .property("class_name", &class_name_wrapper<tcf_type>)
                            .property("module_name", &module_name_wrapper<tcf_type>)
                            .property("sample_name", &sample_name_wrapper<tcf_type>)
                    ]
                ]
            ]
        ]
    ];
}

} // namespace dynamics
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_CORRELATION_HPP */
