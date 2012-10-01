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

#ifndef HALMD_OBSERVABLES_HOST_SAMPLES_BINNED_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_HOST_SAMPLES_BINNED_PHASE_SPACE_HPP

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/clock.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>

namespace halmd {
namespace observables {
namespace host {
namespace samples {

/**
 * store binning data for a phase space sample
 */
template <int dimension, typename float_type>
class binned_phase_space
{
public:
    typedef host::samples::phase_space<dimension, float_type> data_sample_type;
    typedef typename data_sample_type::vector_type vector_type;
    typedef mdsim::clock::step_type step_type;

    /** cell list of particle indices belonging to a single cell/bin */
    typedef std::vector<unsigned int> cell_type;
    /** array of all index cells */
    typedef boost::multi_array<cell_type, dimension> cell_array_type;
    typedef fixed_vector<size_t, dimension> cell_size_type;
    typedef fixed_vector<ssize_t, dimension> cell_diff_type;

    /** d-dimensional array of index cells */
    cell_array_type cell_array;
    /** number of bins per axis, equals shape of cell_array */
    // cell_array.shape() yields 'size_type const*', which is unhandy
    cell_size_type nbin;
    /** edge lengths of a binning cell */
    vector_type cell_length;
    /** origin of lower, left cell (0, … 0) */
    vector_type cell_origin;
    /** simulation step when binning was done */
    step_type step;

    static void luaopen(lua_State* L);

    /**
     * construct binned phase space sample of given size
     *
     * @param data_sample phase space sample underlying the binning
     * @param species particle species selected from sample
     * @param nbin number of bins per space dimension
     */
    binned_phase_space(
        boost::shared_ptr<data_sample_type const> data_sample
      , unsigned int species
      , fixed_vector<unsigned int, dimension> const& nbin
    );

    typename data_sample_type::sample_vector const& position_data() const
    {
        return data_sample_->position(species_);
    }

    typename data_sample_type::sample_vector const& velocity_data() const
    {
        return data_sample_->velocity(species_);
    }

private:
    /** corresponding phase space sample */
    boost::shared_ptr<data_sample_type const> data_sample_;
    /** particle species selected from phase space sample */
    unsigned int species_;
};

} // namespace observables
} // namespace host
} // namespace samples
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_HOST_SAMPLES_BINNED_PHASE_SPACE_HPP */
