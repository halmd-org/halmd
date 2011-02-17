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

#include <boost/array.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/observable.hpp>
#include <halmd/observables/density_modes.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace observables
{

/**
 * compute partial static structure factors
 *
 * @f$ S_q^{(\alpha\beta)} = \langle \frac{1}{N} \rho_{\vec q}^{(\alpha)} \rho_{-\vec q}^{(\beta)} \rangle @f$
 * with the partial density modes
 * @f$ rho_{\vec q}^(\alpha) = \sum_{i=1}^{N_\alpha} \exp(i\vec q\cdot \vec r_i^{\alpha}) @f$
 * and the total number of particles @f$ N = \sum_\alpha N_\alpha @f$
 *
 * see e.g., Hansen & McDonald: Theory of simple liquids, chapter 4.1.
 */

template <int dimension>
class ssf
  : public observable<dimension>
{
public:
    typedef observable<dimension> _Base;
    typedef observables::density_modes<dimension> density_modes_type;
    typedef typename _Base::writer_type writer_type;
    typedef halmd::utility::profiler profiler_type;

    typedef boost::array<double, 3> result_type;
    typedef fixed_vector<double, dimension> vector_type;

    boost::shared_ptr<density_modes_type> density_modes;

    static void luaopen(lua_State* L);

    ssf(
        boost::shared_ptr<density_modes_type> density_modes
      , unsigned int npart
    );

    virtual void register_runtimes(profiler_type& profiler);

    virtual void register_observables(writer_type& writer);

    virtual void prepare() {};

    // compute ssf from sample of density Fourier modes and store with given time
    virtual void sample(double time);

    //! returns last computed values for static structure factor
    std::vector<std::vector<result_type> > const& value() const
    {
        return value_;
    }

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG( sample_, "computation of static structure factor" );

protected:
    /** compute static structure factor and update result accumulators */
    virtual void compute_();

    /** total number of particles, required for normalisation */
    unsigned int npart_;

    /**
     *  result for (partial) static structure factors
     *  in the order AA, AB, AC, …, BB, BC, …,CC
     *
     *  value_[i][j][0]:   mean value S_ab(k) for k = wavenumbers_[j]
     *  value_[i][j][1]:   estimated error of mean
     *  value_[i][j][2]:   value count for the average
     */
    std::vector<std::vector<result_type> > value_;
    /** result accumulators */
    std::vector<std::vector<accumulator<double> > > result_accumulator_;
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
