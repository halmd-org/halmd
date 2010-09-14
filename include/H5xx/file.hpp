/* HDF5 C++ extensions
 *
 * Copyright © 2010  Felix Höfling
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

#ifndef HALMD_UTIL_H5XX_FILE_HPP
#define HALMD_UTIL_H5XX_FILE_HPP

#define H5E_auto_t_vers 2
#include <H5Cpp.h>

namespace H5xx
{

/**
 * HDF5 file
 */
class file : public H5::H5File
{
public:
    /** modes for file creation by constructor or open() */
    enum openmode
    {
        // call H5check() explicitly in the constructors, not via passing a constant!
#       undef H5CHECK
#       define H5CHECK
        // Truncate file, if it already exists, erasing all data previously stored in the file.
        trunc = H5F_ACC_TRUNC,
        // Fail if file already exists. H5F_ACC_TRUNC and H5F_ACC_EXCL are mutually exclusive
        excl = H5F_ACC_EXCL,
        // Open with read only access.
        rdonly = H5F_ACC_RDONLY,
        // Open with read/write access. If the file is currently open for read-only access
        // then it will be reopened.
        rdwr = H5F_ACC_RDWR
    };

    /** default constructor */
    file() { H5check(); }
    /** promote H5::HF5File to H5xx::file */
    file(H5::H5File const& h5file) : H5::H5File(h5file) { H5check(); }

    /**
     * create file in given mode
     */
    file(std::string const& name, openmode mode)
        : H5::H5File(name.c_str(), mode) { H5check(); }
    file(char const* name, openmode mode)
        : H5::H5File(name, mode) { H5check(); }

    /**
     * open file in given mode, only openmode::rdonly or openmode::rdwr are valid
     */
    void open(std::string const& name, openmode mode=rdonly)
    {
        openFile(name.c_str(), mode);
    }

    void open(char const* name, openmode mode=rdonly)
    {
        openFile(name, mode);
    }

};

} // namespace H5xx

#endif /* ! HALMD_UTIL_H5XX_FILE_HPP */
