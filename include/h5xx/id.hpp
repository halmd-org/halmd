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

#include <algorithm>
#include <boost/operators.hpp>
#include <iostream>

#include <h5xx/error.hpp>

namespace h5xx
{

/**
 * This class wraps a HDF5 identifier, which identifies a HDF5
 * library resource, e.g. a file or a group, for safe use within
 * higher-level wrapper classes.
 */
class basic_id
  : boost::equality_comparable<basic_id>
{
public:
    /**
     * wrap identifier returned by H5 API function
     */
    explicit basic_id(hid_t hid)
      : hid_(hid)
    {
    }

    /**
     * copy identifier by increasing reference count
     */
    basic_id(basic_id const& id)
      : hid_(id.hid_)
    {
        H5XX_CHECK(H5Iinc_ref(hid_));
    }

    /**
     * close identifier by decreasing reference count
     */
    ~basic_id()
    {
        H5XX_CHECK(H5Idec_ref(hid_));
    }

    /**
     * assign identifier
     */
    basic_id& operator=(basic_id const& id)
    {
        if (*this != id) {
            basic_id id_(id);
            swap(id_);
        }
        return *this;
    }

    /**
     * returns HDF5 resource identifier
     */
    hid_t hid() const
    {
        return hid_;
    }

    /**
     * swap HDF5 resource identifiers
     */
    void swap(basic_id& id)
    {
        std::swap(hid_, id.hid_);
    }

    /**
     * equality comparison
     */
    bool operator==(basic_id const& id) const
    {
        return hid_ == id.hid_;
    }

    /**
     * returns absolute path of object within file
     *
     * For attributes the path of the parent object is returned.
     */
    std::string path() const
    {
        ssize_t size; // excludes NULL terminator
        H5XX_CHECK(size = H5Iget_name(hid_, NULL, 0));
        std::vector<char> name_(size + 1); // includes NULL terminator
        H5XX_CHECK(H5Iget_name(hid_, name_.data(), name_.size()));
        return name_.data();
    }

private:
    hid_t hid_;
};

/**
 * output absolute path of object within file
 */
std::ostream& operator<<(std::ostream& os, basic_id const& id)
{
    return (os << id.path());
}

} // namespace h5xx

#endif /* ! H5XX_ID_HPP */
