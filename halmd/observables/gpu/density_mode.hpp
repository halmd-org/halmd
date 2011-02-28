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

#ifndef HALMD_OBSERVABLES_GPU_DENSITY_MODE_HPP
#define HALMD_OBSERVABLES_GPU_DENSITY_MODE_HPP

#include <lua.hpp>

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/density_mode.hpp>
#include <halmd/observables/gpu/trajectory.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{
namespace observables { namespace gpu
{

/**
 *  compute Fourier modes of the particle density
 *
 *  @f$ \rho_{\vec q} = \sum_{i=1}^N \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 *  for each particle type
 */
template <int dimension, typename float_type>
class density_mode
  : public observables::density_mode<dimension>
{
public:
    typedef observables::density_mode<dimension> _Base;
    typedef typename _Base::density_mode_sample_type density_mode_sample_type;
    typedef typename _Base::wavevector_type wavevector_type;
    // FIXME use trajectory sample on GPU, add kernels for computation of sum(exp(i q r))
    typedef gpu::trajectory<host::samples::trajectory<dimension, float_type> > trajectory_type;
    typedef halmd::utility::profiler profiler_type;

    typedef fixed_vector<float_type, dimension> vector_type;
    typedef typename density_mode_sample_type::mode_type mode_type;

    static void luaopen(lua_State* L);

    density_mode(
        boost::shared_ptr<trajectory_type> trajectory
      , boost::shared_ptr<wavevector_type> wavevector
    );

    void register_runtimes(profiler_type& profiler);

    /**
    * compute density modes from trajectory sample and store with given time stamp
    */
    virtual void acquire(double time);

    //! returns nested list of density modes
    virtual typename _Base::result_type const& value() const
    {
        return rho_sample_.rho;
    }

    //! returns wavevector object
    virtual wavevector_type const& wavevector() const
    {
        return *wavevector_;
    }

    //! returns wavenumber grid
    virtual std::vector<double> const& wavenumber() const
    {
        return wavevector_->wavenumber();
    }

    // descriptions of module's runtime accumulators
    HALMD_PROFILING_TAG(sample_, "computation of density modes");

protected:
    boost::shared_ptr<trajectory_type> trajectory_;
    boost::shared_ptr<wavevector_type> wavevector_;

    /** data structure for density modes */
    density_mode_sample_type rho_sample_;

    // list of profiling timers
    boost::fusion::map<
        boost::fusion::pair<sample_, accumulator<double> >
    > runtime_;
};

}} // namespace observables::gpu

}  // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DENSITY_MODE_HPP */
