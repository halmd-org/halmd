/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/io/trajectory/reader.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory
{

/**
 * Assemble module options
 */
template <int dimension>
void reader<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("trajectory-file,J", po::value<string>()->required(),
         "trajectory input file")
        ;

    po::options_description group("Trajectory reader");
    group.add_options()
        ("trajectory-sample,S", po::value<ssize_t>()->required(),
         "trajectory sample for initial state")
        ;
    desc.add(group);
}

template <int dimension>
reader<dimension>::reader(po::options const& vm)
  // parse options
  : path_(vm["trajectory-file"].as<string>())
  , offset_(vm["trajectory-sample"].as<ssize_t>())
{
}

template class reader<3>;
template class reader<2>;

}} // namespace io::trajectory

} // namespace halmd
