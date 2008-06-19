/* Time statistics
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

#ifndef MDSIM_PERFORMANCE_HPP
#define MDSIM_PERFORMANCE_HPP

#include <H5Cpp.h>
#include <iomanip>
#include <map>
#include <string>
#include "H5param.hpp"
#include "accumulator.hpp"
#include "log.hpp"
#include "statistics.hpp"


namespace mdsim
{

/**
 * time accumulators
 */
typedef std::map<std::string, accumulator<double> > perf_type;

/**
 * time statistics
 */
template <unsigned dimension, typename T>
class perf
{
public:
    perf();

    /** create HDF5 time statistics output file */
    void open(std::string const& filename);
    /** dump global simulation parameters to HDF5 file */
    perf<dimension, T>& operator<<(H5param const& param);
    /** write time statistics to HDF5 file */
    void write(perf_type const& times);
    /** close HDF5 file */
    void close();

private:
    /** HDF5 time statistics output file */
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
 * create HDF5 time statistics output file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::open(std::string const& filename)
{
    LOG("write time statistics to file: " << filename);
    try {
	// truncate existing file
	file_ = H5::H5File(filename, H5F_ACC_TRUNC);
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to create time statistics file");
    }
}

/**
 * dump global simulation parameters to HDF5 file
 */
template <unsigned dimension, typename T>
perf<dimension, T>& perf<dimension, T>::operator<<(H5param const& param)
{
    param.write(file_.createGroup("/parameters"));
    return *this;
}

/**
 * write time statistics to HDF5 file
 */
template <unsigned dimension, typename T>
void perf<dimension, T>::write(perf_type const& times)
{
    hsize_t dim[1] = { 3 };
    H5::DataSpace ds(1, dim);
    H5::DataType tid(H5::PredType::NATIVE_DOUBLE);
    H5::Group root(file_.createGroup("/times"));

    try {
	for (typename perf_type::const_iterator it = times.begin(); it != times.end(); ++it) {
	    // create dataset
	    H5::DataSet dataset(root.createDataSet(it->first, tid, ds));

	    // write dataset
	    double data[] = {
		// average time in milliseconds
		it->second.mean(),
		// standard deviation in milliseconds
		it->second.std(),
		// number of calls
		it->second.count(),
	    };
	    dataset.write(data, tid, ds, ds);

	    if (it->second.count() > 1) {
		LOG("average " << it->first << " time: " << std::setprecision(4) << data[0] << " ms (" << std::setprecision(4) << data[1] << " ms)");
	    }
	    else {
		LOG(it->first << " time: " << std::setprecision(4) << data[0] << " ms");
	    }
	}
    }
    catch (H5::FileIException const& e) {
	throw exception("failed to write time statistics to HDF5 file");
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
	throw exception("failed to close time statistics file");
    }
}

} // namespace mdsim

#endif /* ! MDSIM_PERFORMANCE_HPP */
