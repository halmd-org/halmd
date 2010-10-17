/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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

#ifndef HALMD_TEST_TOOLS_CUDA_HPP
#define HALMD_TEST_TOOLS_CUDA_HPP

#include <boost/test/unit_test_monitor.hpp>
#include <iostream>

#include <cuda_wrapper/cuda_wrapper.hpp>

// "global fixture:" select CUDA device
struct set_cuda_device {
    set_cuda_device() {
        try {
            cuda::device::set(0);
        }
        catch (cuda::error const& e) {
            std::cerr << "CUDA error: " << e.what() << std::endl;
            throw;
        }
    }
    ~set_cuda_device() {}  // release device here?
};

void cuda_error_translator( cuda::error const& e)
{
    BOOST_TEST_MESSAGE( "(CUDA error) " << e.what() ); throw;
}

struct register_cuda_error {
    register_cuda_error() {
        using namespace boost::unit_test;
        unit_test_monitor.register_exception_translator<cuda::error>( &cuda_error_translator );
    }
};

static register_cuda_error _register_cuda_error_dummy;

#endif /* ! HALMD_TEST_TOOLS_CUDA_HPP */
