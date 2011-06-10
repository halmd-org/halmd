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

#ifndef HALMD_OBSERVABLES_DYNAMICS_TCF_HPP
#define HALMD_OBSERVABLES_DYNAMICS_TCF_HPP

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>

namespace halmd
{
namespace observables { namespace dynamics
{

/**
 * Compute a time correlation function of data stored
 * according to the 'blocking scheme'.
 */
class tcf_base
{
public:
    tcf_base() {}
    virtual ~tcf_base() {}

    /** compute correlations at the given coarse-graining level */
    virtual void compute(unsigned int level) = 0;
};

template <typename Functor>
class tcf
  : public tcf_base
{
public:
    typedef typename Functor::sample_type sample_type;
    typedef typename Functor::result_type result_type;
    typedef observables::samples::blocking_scheme<sample_type> block_sample_type;
    typedef boost::multi_array<accumulator<result_type>, 2> block_result_type;

    static void luaopen(lua_State* L, char const* scope, char const* class_name);

    tcf(
        boost::shared_ptr<Functor> tcf
      , boost::shared_ptr<block_sample_type> block_sample
    );
    virtual ~tcf() {}

    virtual void compute(unsigned int level);

private:
    typedef tcf_base _Base;

    /** block structure holding the input data */
    boost::shared_ptr<block_sample_type> block_sample_;
    /** functor performing the specific computation */
    boost::shared_ptr<Functor> tcf_;

    /** block structures holding accumulated result values */
    block_result_type result_;
};

template <typename Functor>
tcf<Functor>::tcf(
    boost::shared_ptr<Functor> tcf
  , boost::shared_ptr<block_sample_type> block_sample
)
  // dependency injection
  : block_sample_(block_sample)
  , tcf_(tcf)
  // memory allocation
  , result_(boost::extents[block_sample->count()][block_sample->block_size()])
{
}

template <typename Functor>
void tcf<Functor>::compute(unsigned int level)
{
    LOG_TRACE("[" << Functor::module_name() << "]: compute correlations at level " << level);

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

template <typename Functor>
static char const* module_name_wrapper(tcf<Functor> const&)
{
    return Functor::module_name();
}

template <typename Functor>
static char const* class_name_wrapper(tcf<Functor> const&)
{
    return Functor::class_name();
}

template <typename Functor>
void tcf<Functor>::luaopen(lua_State* L, char const* scope, char const* class_name)
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
                    namespace_("tcf")
                    [
                        class_<tcf, boost::shared_ptr<_Base>, _Base>(class_name)
                            .def(constructor<
                                boost::shared_ptr<Functor>
                              , boost::shared_ptr<block_sample_type>
                            >())
                            .property("module_name", &module_name_wrapper<Functor>)
                            .property("class_name", &class_name_wrapper<Functor>)
                    ]
                ]
            ]
        ]
    ];
}

}} // namespace observables::dynamics

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_DYNAMICS_TCF_HPP */
