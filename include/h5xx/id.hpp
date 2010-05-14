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

#ifndef H5XX_ID_HPP
#define H5XX_ID_HPP

#include <boost/operators.hpp>

#include <h5xx/error.hpp>

namespace h5xx
{

/**
 * HDF5 identifier
 */
class id
  : public boost::equality_comparable<id>
{
public:
    /**
     * close identifier
     */
    ~id()
    {
        H5XX_CALL(H5Idec_ref(id_));
    }

    /**
     * copy identifier
     */
    id(id const& other)
    {
        H5XX_CALL(H5Iinc_ref(id_ = other.id_));
    }

    /**
     * assign identifier
     */
    id& operator=(id const& other)
    {
        if (id_ != other.id_) {
            H5XX_CALL(H5Idec_ref(id_));
            H5XX_CALL(H5Iinc_ref(id_ = other.id_));
        }
        return *this;
    }

    /**
     * equality comparison
     */
    bool operator==(id const& other) const
    {
        return id_ == other.id_;
    }

protected:
    explicit id(hid_t id_) : id_(id_) {}
    hid_t id_;
};

} // namespace h5xx

#endif /* ! H5XX_ID_HPP */
