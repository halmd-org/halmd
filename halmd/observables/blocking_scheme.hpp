/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_BLOCKING_SCHEME_HPP
#define HALMD_OBSERVABLES_BLOCKING_SCHEME_HPP

#include <boost/shared_ptr.hpp>

#include <halmd/observables/samples/blocking_scheme.hpp>

namespace halmd
{
namespace observables
{

class blocking_scheme
{
public:
    typedef samples::blocking_scheme_base block_sample_type;

    blocking_scheme(boost::shared_ptr<block_sample_type> block_sample);

private:
    boost::shared_ptr<block_sample_type> block_sample_;
};

} // namespace observables

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_BLOCKING_SCHEME_HPP */
