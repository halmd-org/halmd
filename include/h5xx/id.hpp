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
#include <iostream>

#include <h5xx/error.hpp>

namespace h5xx
{

/**
 * HDF5 identifier
 */
class id
  : boost::equality_comparable<id>
{
public:
    /**
     * close identifier
     */
    ~id()
    {
        H5XX_CHECK(H5Idec_ref(id_));
    }

    /**
     * copy identifier
     */
    id(id const& other)
    {
        H5XX_CHECK(H5Iinc_ref(id_ = other.id_));
    }

    /**
     * assign identifier
     */
    id& operator=(id const& other)
    {
        if (id_ != other.id_) {
            H5XX_CHECK(H5Idec_ref(id_));
            H5XX_CHECK(H5Iinc_ref(id_ = other.id_));
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

    /**
     * returns absolute path of object within file
     */
    std::string path() const
    {
        ssize_t size; // excludes NULL terminator
        H5XX_CHECK(size = H5Iget_name(id_, NULL, 0));
        std::vector<char> name_(size + 1); // includes NULL terminator
        H5XX_CHECK(H5Iget_name(id_, name_.data(), name_.size()));
        return name_.data();
    }

protected:
    explicit id(hid_t id_)
      : id_(id_)
    {}

    hid_t id_;
};

/**
 * output absolute path of object within file
 */
std::ostream& operator<<(std::ostream& os, id const& id_)
{
    return (os << id_.path());
}

} // namespace h5xx

#endif /* ! H5XX_ID_HPP */
