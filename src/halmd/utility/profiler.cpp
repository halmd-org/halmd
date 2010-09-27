/*
 * Copyright © 2010  Peter Colberg
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

#include <halmd/io/profile/writer.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/profiler.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace utility
{

/**
 * Resolve module dependencies
 */
void profiler::depends()
{
    modules::depends<_Self, profile_writer_type>::required();
}

profiler::profiler(modules::factory& factory, po::variables_map const& vm)
  // dependency injection
  : profile_writer(modules::fetch<profile_writer_type>(factory, vm)) {}

void profiler::register_accumulator(
    std::type_info const& tag
  , accumulator_type const& acc
  , std::string const& desc
) const
{
    for_each(
        profile_writer.begin()
      , profile_writer.end()
      , bind(
            &profile_writer_type::register_accumulator
          , _1
          , tokenized_name(tag) // namespace and class template tokens
          , cref(acc)
          , cref(desc)
        )
    );
}

} // namespace utility

template class module<utility::profiler>;

} // namespace halmd
