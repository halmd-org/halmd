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

#include <fstream>

#include <halmd/mdsim/random.hpp>
#include <halmd/util/log.hpp>

using namespace std;

namespace halmd { namespace mdsim
{

/**
 * Obtain integer seed from file
 */
unsigned int random::readint(std::string const& file)
{
    LOG("reading 32-bit integer seed from " << file);
    unsigned int seed;
    try {
        ifstream f;
        f.exceptions(ifstream::eofbit|ifstream::failbit|ifstream::badbit);
        f.open(file.c_str());
        f.read(reinterpret_cast<char*>(&seed), sizeof(seed));
        f.close();
    }
    catch (ifstream::failure const& e) {
        throw logic_error("failed to read from " + file);
    }
    return seed;
}

}} // namespace halmd::mdsim
