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

#ifndef H5XX_FILE_HPP
#define H5XX_FILE_HPP

#include <boost/preprocessor/tuple/elem.hpp>
#include <string>
#include <vector>

#include <h5xx/group.hpp>
#include <h5xx/id.hpp>

namespace h5xx
{

/**
 * HDF5 file
 */
class file
  : public virtual basic_id
  , public basic_group
{
public:
// The H5 file mode macros cannot be used to define an enum,
// as these macros call H5check() as a side effect. We apply
// some boost preprocessor magic to extract the flag.

#define H5XX_FILE_FLAG(x) BOOST_PP_TUPLE_ELEM(2, 1, x)

    enum open_mode {
        /** open in read-only mode */
        ro = H5XX_FILE_FLAG(H5F_ACC_RDONLY),
        /** open in read-write mode */
        rw = H5XX_FILE_FLAG(H5F_ACC_RDWR),
    };

    enum create_mode {
        /** truncate existing file */
        trunc = H5XX_FILE_FLAG(H5F_ACC_TRUNC),
        /** fail if file does not exist */
        excl  = H5XX_FILE_FLAG(H5F_ACC_EXCL),
    };

    /**
     * open existing HDF5 file
     */
    file(std::string const& name, open_mode flags = file::ro)
      : basic_id(open(name, flags))
      , basic_group(basic_id::hid())
    {}

    /**
     * create and open HDF5 file
     */
    file(std::string const& name, create_mode flags)
      : basic_id(create(name, flags))
      , basic_group(basic_id::hid())
    {}

public:
    /**
     * flushes all buffers associated with file to disk
     */
    void flush()
    {
        H5XX_CHECK(H5Fflush(basic_id::hid(), H5F_SCOPE_LOCAL));
    }

    /**
     * returns file name
     */
    std::string name()
    {
        ssize_t size; // excludes NULL terminator
        H5XX_CHECK(size = H5Fget_name(basic_id::hid(), NULL, 0));
        std::vector<char> name_(size + 1); // includes NULL terminator
        H5XX_CHECK(H5Fget_name(basic_id::hid(), name_.data(), name_.size()));
        return name_.data();
    }

    /**
     * checks whether a file is in the HDF5 format
     */
    static bool format(std::string const& name)
    {
        htri_t result; // no bool, otherwise error handling will fail
        H5XX_CHECK(result = H5Fis_hdf5(name.c_str()));
        return result;
    }

private:
    /**
     * open existing HDF5 file
     */
    static basic_id open(std::string const& name, unsigned int flags)
    {
        H5XX_CHECK(H5check()); // compare header and library versions
        hid_t hid;
        H5XX_CHECK(hid = H5Fopen(name.c_str(), flags, H5P_DEFAULT));
        return basic_id(hid);
    }

    /**
     * create and open HDF5 file
     */
    static basic_id create(std::string const& name, unsigned int flags)
    {
        H5XX_CHECK(H5check());
        hid_t hid;
        H5XX_CHECK(hid = H5Fcreate(name.c_str(), flags, H5P_DEFAULT, H5P_DEFAULT));
        return basic_id(hid);
    }
};

} // namespace h5xx

#endif /* ! H5XX_FILE_HPP */
