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

#include <boost/assign.hpp>
#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <string>
#include <ljgpu/sample/H5param.hpp>
#include <ljgpu/math/accum.hpp>
#include <ljgpu/util/H5xx.hpp>

namespace ljgpu
{

/**
 * performance data
 */
class perf : boost::noncopyable
{
public:
    /* performance accumulators */
    typedef boost::unordered_map<std::string, accumulator<float> > counters;
    /* performance accumulator */
    typedef counters::value_type counter;
    /* performance counter descriptions */
    typedef boost::unordered_map<std::string, std::string> desc_map;

public:
    perf() : m_is_open(false) {}
    /** create HDF5 performance data output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    operator H5param() { return m_file; }
    /** sample performance data */
    void sample(counters const& times);
    /** write performance data to HDF5 file */
    void flush();
    /** close HDF5 file */
    void close();
    /** returns true iff associated with HDF5 file */
    bool is_open() const { return m_is_open; }

    /** returns accumulators */
    counters const& times() const { return m_times; }
    /** returns accumulator description */
    static std::string const& desc(std::string name) { return m_desc.at(name); }

private:
    /** CPU tick accumulators */
    counters m_times;
    /** HDF5 performance data output file */
    H5::H5File m_file;
    /** true iff associated with HDF5 file */
    bool m_is_open;
    /** HDF5 datasets */
    boost::unordered_map<std::string, H5::DataSet> m_dataset;
    /** HDF5 floating-point data type */
    H5::DataType m_tid;
    /* performance counter descriptions */
    static desc_map m_desc;
};

/** output formatted performance statistics to stream */
std::ostream& operator<<(std::ostream& os, perf::counters const& times);

} // namespace ljgpu

#endif /* ! LJGPU_SAMPLE_PERF_HPP */
