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

#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include <set>
#include <sstream>

#include <halmd/utility/module/options.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace utility { namespace module
{

/**
 * outputs options description of all modules
 */
ostream& operator<<(ostream& os, options_description const& opt)
{
    // First we collect the options of all modules and output
    // them to a string stream. Options not within a separate
    // options group are listed in the Modules section.

    po::options_description desc("Modules");
    BOOST_FOREACH( factory::_Module_map::value_type const& pair, factory::modules() ) {
        pair.second->options(desc);
    }
    stringstream s;
    s << desc;

    // Group options by sections and eliminate duplicate options,
    // using double-ended queues to maintain insertion order.

    pair<deque<string>, map<string, pair<deque<string>, set<string> > > > opts;
    string line, section;
    while (getline(s, line)) {
        if (ends_with(line, ":")) { // start of new section
            section = line;
            if (!opts.second.count(section)) {
               opts.first.push_back(section);
            }
        }
        else if (!line.empty()) {
            if (opts.second[section].second.insert(line).second) {
                opts.second[section].first.push_back(line);
            }
        }
    }

    // Output unique option groups.

    BOOST_FOREACH( string const& section, opts.first ) {
        os << section << endl;
        BOOST_FOREACH( string const& option, opts.second[section].first ) {
            os << option << endl;
        }
        os << endl;
    }
    return os;
}

}} // namespace utility::module

} // namespace halmd
