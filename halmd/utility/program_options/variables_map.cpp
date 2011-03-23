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

#include <boost/foreach.hpp>
#include <halmd/utility/program_options/errors.hpp>
#include <halmd/utility/program_options/typed_value.hpp>
#include <halmd/utility/program_options/variables_map.hpp>

using namespace boost;
using namespace boost::program_options;
using namespace std;

namespace halmd
{
namespace po
{

/**
 * extends boost::program_options::store with conflicting/dependent options
 */
void store(parsed_options const& options, variables_map& vm, bool utf8)
{
    program_options::store(options, vm, utf8);

    options_description const& desc = *options.description;

    BOOST_FOREACH(shared_ptr<option_description> const& d, desc.options())
    {
        string const& name = d->long_name();

        // cast upstream value semantic to our value semantic
        shared_ptr<extended_value_semantic const> semantic(
            dynamic_pointer_cast<extended_value_semantic const>(d->semantic())
        );
        // skip values not deriving from our value semantic
        if (!semantic) {
            continue;
        }

        // validate conflicting/dependent options exist
        BOOST_FOREACH( string const& second, semantic->conflicting() ) {
#if BOOST_VERSION >= 104200
            desc.find(second, false, false, false); //< throws unknown_option
#else
            desc.find(second, false);
#endif
        }
        BOOST_FOREACH( string const& second, semantic->dependent() ) {
#if BOOST_VERSION >= 104200
            desc.find(second, false, false, false);
#else
            desc.find(second, false);
#endif
        }

        // skip not parsed or defaulted option
        if (vm[name].empty() || vm[name].defaulted()) {
            continue;
        }

        // validate option conflicts
        BOOST_FOREACH( string const& other, semantic->conflicting() ) {
            // throw if parsed and non-defaulted option
            if (!vm[other].empty() && !vm[other].defaulted()) {
                throw conflicting_option(name, other);
            }
        }
        // validate option dependencies
        BOOST_FOREACH( string const& other, semantic->dependent() ) {
            // throw if not parsed or defaulted option
            if (vm[other].empty() || vm[other].defaulted()) {
                throw dependent_option(name, other);
            }
        }
    }
}

void store(wparsed_options const& options, variables_map& vm)
{
    po::store(options.utf8_encoded_options, vm, true);
}

} // namespace po

} // namespace halmd
