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

#ifndef HALMD_OBSERVABLES_HOST_DENSITY_MODE_HPP
#define HALMD_OBSERVABLES_HOST_DENSITY_MODE_HPP

#include <lua.hpp>

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/density_mode.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace host {

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
    typedef host::samples::phase_space<dimension, float_type> phase_space_type;
    typedef halmd::utility::profiler profiler_type;
    typedef typename _Base::signal_type signal_type;
    typedef typename _Base::slot_function_type slot_function_type;
    typedef typename _Base::connection_type connection_type;

    typedef fixed_vector<float_type, dimension> vector_type;
    typedef typename density_mode_sample_type::mode_type mode_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type sample;
    };

    static void luaopen(lua_State* L);

    density_mode(
        boost::shared_ptr<phase_space_type const> phase_space
      , boost::shared_ptr<wavevector_type const> wavevector
    );

    void register_runtimes(profiler_type& profiler);

    /**
    * compute density modes from phase space sample and store with given time stamp
    */
    virtual void acquire(uint64_t step);

    virtual connection_type on_acquire(slot_function_type const& slot)
    {
        return on_acquire_.connect(slot);
    }

    //! returns nested list of density modes
    virtual typename _Base::result_type const& value() const
    {
        return rho_sample_.rho;
    }

    //! returns simulation step when sample was taken
    virtual uint64_t step() const
    {
        return rho_sample_.step;
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

private:
    boost::shared_ptr<phase_space_type const> phase_space_;
    boost::shared_ptr<wavevector_type const> wavevector_;

    /** data structure for density modes */
    density_mode_sample_type rho_sample_;
    /** profiling runtime accumulators */
    runtime runtime_;

    signal_type on_acquire_;
};

} // namespace observables
} // namespace host
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_DENSITY_MODE_HPP */
