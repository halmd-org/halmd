/* Performance data
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_SAMPLE_PERF_HPP
#define LJGPU_SAMPLE_PERF_HPP

#include <H5Cpp.h>
#include <boost/assign.hpp>
#include <boost/unordered_map.hpp>
#include <string>
#include <ljgpu/util/H5param.hpp>
#include <ljgpu/math/accum.hpp>

namespace ljgpu
{

/**
 * performance data
 */
class perf
{
public:
    /* performance accumulators */
    typedef boost::unordered_map<std::string, accumulator<float> > counters;
    /* performance accumulator */
    typedef counters::value_type counter;
    /* performance counter descriptions */
    typedef boost::unordered_map<std::string, std::string> desc_map;

public:
    perf() : m_offset(0), m_dirty(false) {}
    /** create HDF5 performance data output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** sample performance data */
    void sample(counters const& times);
    /** clear performance counters */
    void commit();
    /** write performance data to HDF5 file */
    void flush(bool force=true);
    /** close HDF5 file */
    void close();

private:
    void create_datasets();
    void write_datasets();

private:
    /** CPU tick accumulators */
    counters m_times;
    /** HDF5 performance data output file */
    H5::H5File m_file;
    /** HDF5 datasets */
    boost::unordered_map<std::string, H5::DataSet> m_dataset;
    /** HDF5 floating-point data type */
    H5::DataType m_tid;
    /** dataset offset */
    uint64_t m_offset;
    /** pending data bit */
    bool m_dirty;
    /* performance counter descriptions */
    static desc_map desc;
};

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_PERF_HPP */
