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
#include <halmd/utility/lua.hpp>

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
        ("trajectory-file,J", po::value<string>(),
         "trajectory input file")
        ;

    po::options_description group("Trajectory reader");
    group.add_options()
        ("trajectory-sample,S", po::value<ssize_t>()->default_value(-1),
         "trajectory sample for initial state")
        ;
    desc.add(group);
}

template <int dimension>
void reader<dimension>::select(po::variables_map const& vm)
{
    if (vm["trajectory-file"].empty()) {
        throw unsuitable_module("mismatching option trajectory-file");
    }
}

template <int dimension>
reader<dimension>::reader(modules::factory& factory, po::variables_map const& vm)
  // parse options
  : path_(vm["trajectory-file"].as<string>())
  , offset_(vm["trajectory-sample"].as<ssize_t>())
{
}

template <typename T>
static luabind::scope register_lua(char const* class_name)
{
    using namespace luabind;
    return
        namespace_("halmd_wrapper")
        [
            namespace_("io")
            [
                namespace_("trajectory")
                [
                    class_<T, shared_ptr<T> >(class_name)
                        .scope
                        [
                            def("options", &T::options)
                        ]
                ]
            ]
        ];
}

static lua_registry::iterator dummy = (
    lua_registry::get()->push_back( register_lua<reader<3> >("reader_3_") )
  , lua_registry::get()->push_back( register_lua<reader<2> >("reader_2_") )
  , lua_registry::get()->end()
);

template class reader<3>;
template class reader<2>;

}} // namespace io::trajectory

} // namespace halmd
