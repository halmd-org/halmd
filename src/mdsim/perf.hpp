/* Performance data
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_PERF_HPP
#define MDSIM_PERF_HPP

#include <H5Cpp.h>
#include <boost/array.hpp>
#include <string>
#include "H5param.hpp"
#include "accumulator.hpp"

namespace mdsim
{

/**
 * performance class accumulators
 */
typedef boost::array<accumulator<double>, 5> perf_counters;

/**
 * performance data
 */
class perf
{
public:
    perf() : m_offset(0), m_dirty(false) {}
    /** create HDF5 performance data output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** sample performance data */
    void sample(perf_counters const& times);
    /** clear performance counters */
    void commit();
    /** write performance data to HDF5 file */
    void flush(bool force=true);
    /** close HDF5 file */
    void close();

private:
    /** CPU tick accumulators */
    perf_counters m_times;
    /** HDF5 performance data output file */
    H5::H5File m_file;
    /** HDF5 datasets */
    boost::array<H5::DataSet, 5> m_dataset;
    /** HDF5 floating-point data type */
    H5::DataType m_tid;
    /** dataset offset */
    uint64_t m_offset;
    /** pending data bit */
    bool m_dirty;
};

} // namespace mdsim

#endif /* ! MDSIM_PERF_HPP */
