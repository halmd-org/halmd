/* MD simulation trajectory writer
 *
 * Copyright © 2008-2009  Peter Colberg
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

#include <halmd/sample/trajectory.hpp>

namespace halmd {

/**
 * create HDF5 trajectory output file
 */
void trajectory::open(std::string const& filename, openmode mode)
{
    if (mode & out) {
        LOG("write trajectories to file: " << filename);
        try {
            // truncate existing file
            m_file = H5::H5File(filename, H5F_ACC_TRUNC);
            m_is_open = true;
        }
        catch (H5::FileIException const& e) {
            throw exception("failed to create trajectory output file");
        }
    }
    else {
        LOG("read trajectory file: " << filename);
        try {
            m_file = H5::H5File(filename, H5F_ACC_RDONLY);
            m_is_open = true;
        }
        catch (H5::Exception const& e) {
            throw exception("failed to open HDF5 trajectory input file");
        }
    }
}

/**
 * close HDF5 trajectory file
 */
void trajectory::close()
{
    try {
        m_file.close();
        m_is_open = false;
    }
    catch (H5::Exception const& e) {
        throw exception("failed to close HDF5 trajectory file");
    }
}

/**
 * flush HDF5 output file to disk
 */
void trajectory::flush()
{
    try {
        m_file.flush(H5F_SCOPE_GLOBAL);
    }
    catch (H5::Exception const& e) {
        throw exception("failed to flush HDF5 trajectory output file");
    }
}

} // namespace halmd
