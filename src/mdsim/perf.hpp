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
#include <iomanip>
#include <map>
#include <string>
#include <sys/times.h>
#include <unistd.h>
#include "H5param.hpp"
#include "accumulator.hpp"
#include "log.hpp"
#include "statistics.hpp"


namespace mdsim
{

/**
 * performance class accumulators
 */
typedef std::map<std::string, std::map<std::string, accumulator<double> > > perf_type;

/**
 * performance data
 */
template <unsigned dimension, typename T>
class perf
{
public:
    perf();

    /** create HDF5 performance data output file */
    void open(std::string const& filename);
    /** returns HDF5 parameter group */
    H5param attrs();
    /** write performance data to HDF5 file */
    void write(perf_type const& times);
    /** close HDF5 file */
    void close();

private:
    /** HDF5 performance data output file */
    H5::H5File file_;
};


template <unsigned dimension, typename T>
perf<dimension, T>::perf()
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif
}

/**
 * create HDF5 performance data output file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::open(std::string const& filename)
{
    LOG("write performance data to file: " << filename);
    try {
	// truncate existing file
	file_ = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create performance data file");
    }
    // create parameter group
    file_.createGroup("param");
}

/**
 * returns HDF5 parameter group
 */
template <unsigned dimension, typename T>
H5param perf<dimension, T>::attrs()
{
    return H5param(file_.openGroup("param"));
}

/**
 * write performance data to HDF5 file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::write(perf_type const& times)
{
    H5::DataType ftid(H5::PredType::NATIVE_DOUBLE);
    H5::DataType utid(H5::PredType::NATIVE_UINT64);

    // milliseconds per clock tick
    const double tick = 1. / sysconf(_SC_CLK_TCK);

    try {
	H5::Group root(file_.createGroup("/times"));

	for (typename perf_type::const_iterator it = times.begin(); it != times.end(); ++it) {
	    // create group for performance class
	    H5::Group group(root.createGroup(it->first));

	    for (typename perf_type::mapped_type::const_iterator acc = it->second.begin(); acc != it->second.end(); ++acc) {
		// average time in seconds
		const double mean = acc->second.mean() * tick;
		// standard deviation in seconds
		const double sigma = acc->second.std() * tick;
		// number of calls
		const uint64_t count = acc->second.count();

		H5::Group time(group.createGroup(acc->first));
		time.createAttribute("mean", ftid, H5S_SCALAR).write(ftid, &mean);
		time.createAttribute("sigma", ftid, H5S_SCALAR).write(ftid, &sigma);
		time.createAttribute("count", utid, H5S_SCALAR).write(utid, &count);

		if (acc->second.count() > 1) {
		    LOG(acc->first << "(" << it->first << ") average time: " << std::setprecision(4) << (mean * 1.E3) << " ms (" << std::setprecision(4) << (sigma * 1.E3) << " ms, " << count << " calls)");
		}
		else {
		    LOG(acc->first << "(" << it->first << ") time: " << std::setprecision(4) << (mean * 1.E3) << " ms");
		}
	    }
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write performance data to HDF5 file");
    }
}

/**
 * close HDF5 file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::close()
{
    try {
	file_.close();
    }
    catch (H5::Exception const& e) {
	throw exception("failed to close performance data file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_PERF_HPP */
