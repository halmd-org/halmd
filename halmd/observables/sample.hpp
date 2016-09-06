/*
 * Copyright 2016 Daniel Kirchner
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HALMD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HALMD.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_OBSERVABLES_SAMPLE_HPP
#define HALMD_OBSERVABLES_SAMPLE_HPP

#include <typeinfo>

namespace halmd {
namespace observables {

/** base class for host and GPU samples */
class sample_base {
public:
    virtual ~sample_base() {}
    virtual std::type_info const& type() const = 0;
};

/** helper structure to identify GPU samples
 * sample_base::type() will return typeid(gpu_sample<T>) for a GPU sample of type T
 */
template<typename T>
struct gpu_sample {
    typedef T type;
};

} // namespace observables
} // namespace halmd

#endif /* !defined HALMD_OBSERVABLES_SAMPLE_HPP */
