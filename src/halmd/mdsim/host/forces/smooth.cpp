/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <boost/algorithm/string/predicate.hpp>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;

namespace halmd
{
namespace mdsim { namespace host { namespace forces
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void smooth<dimension, float_type>::options(po::options_description& desc)
{
    desc.add_options()
        ("smooth", po::value<float>()->default_value(0.f),
         "C²-potential smoothing factor")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void smooth<dimension, float_type>::resolve(po::options const& vm)
{
    using namespace boost;
    if (!ends_with(vm["force"].as<std::string>(), "-smooth")) {
        throw unsuitable_module<smooth>("mismatching option '--force'");
    }
}

/**
 * Initialise parameters
 */
template <int dimension, typename float_type>
smooth<dimension, float_type>::smooth(po::options const& vm)
  // initialise parameters
  : r_smooth_(vm["smooth"].as<float>())
  , rri_smooth_(std::pow(r_smooth_, -2))
{
    LOG("scale parameter for potential smoothing function: " << r_smooth_);
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class smooth<3, double>;
template class smooth<2, double>;
#else
template class smooth<3, float>;
template class smooth<2, float>;
#endif

}}} // namespace mdsim::host::forces

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::forces::smooth<3, double> >;
template class module<mdsim::host::forces::smooth<2, double> >;
#else
template class module<mdsim::host::forces::smooth<3, float> >;
template class module<mdsim::host::forces::smooth<2, float> >;
#endif

} // namespace halmd
