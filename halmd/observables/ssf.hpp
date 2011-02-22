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

#ifndef HALMD_OBSERVABLES_SSF_HPP
#define HALMD_OBSERVABLES_SSF_HPP

#include <boost/tuple/tuple.hpp>
#include <lua.hpp>
#include <map>
#include <vector>

#include <halmd/observables/observable.hpp>
#include <halmd/observables/samples/density_modes.hpp>
#include <halmd/observables/utility/wavevectors.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace observables
{

/**
 * compute static structure factor
 */

template <int dimension>
class ssf
  : public observable<dimension>
{
public:
    typedef observable<dimension> _Base;

    typedef samples::density_modes<dimension> density_modes_type;
    typedef utility::wavevectors<dimension> wavevectors_type;
    typedef typename _Base::writer_type writer_type;
    typedef halmd::utility::profiler profiler_type;
    typedef fixed_vector<double, dimension> vector_type;

    boost::shared_ptr<density_modes_type> density_modes;
    boost::shared_ptr<wavevectors_type> wavevectors;

    static void luaopen(lua_State* L);

    ssf(
        boost::shared_ptr<density_modes_type> density_modes
      , boost::shared_ptr<wavevectors_type> wavevectors
    );

    void register_runtimes(profiler_type& profiler);
    virtual void register_observables(writer_type& writer);

    virtual void prepare() {}

    // compute ssf from trajectory sample and store with given time
    virtual void sample(double time);

    //! returns last computed result for static structure factor
    std::vector<boost::tuple<double, double, unsigned> > const& result() const
    {
        return result_;
    }

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG( sample_, "computation of static structure factor" );

protected:
    /** compute static structure factor and update accumulators 'result_acc_' */
    virtual void compute_() = 0;

    /**
     *  result for static structure factor
     *
     *  result_[i][0]:   mean value S(k) for k = wavenumbers_[i]
     *  result_[i][1]:   estimated error of mean
     *  result_[i][2]:   value count for the average
     */
    std::vector<boost::tuple<double, double, unsigned> > result_;
    // FIXME support HDF5 output of tuples
    std::vector<double> value_;
    std::vector<double> error_;
    std::vector<unsigned int> count_;
    /** result accumulators */
    std::vector<accumulator<double> > result_acc_;
    /** store time to be passed to HDF5 writer */
    double time_;

    // list of profiling timers
    boost::fusion::map<
        boost::fusion::pair<sample_, accumulator<double> >
    > runtime_;
};

} // namespace observables

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SSF_HPP */
