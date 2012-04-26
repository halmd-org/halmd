/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test_monitor.hpp>
#include <iostream>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <test/tools/ctest.hpp>

/**
 *  "global fixture" for Boost Unit Test Framework: select CUDA device
 *
 *  create a CUDA context on a free device using the driver library, respect device locks
 */
class set_cuda_device
{
public:
    set_cuda_device() {
        // choose first available CUDA device
        int device_count = cuda::device::count();
        for (int i = 0; i < device_count; ++i) {
            try {
                // create CUDA context and associate it with this thread
                ctx_.reset(new cuda::driver::context(i));
                break;
            }
            catch (cuda::driver::error const&) {
                // device is compute-exlusive mode and in use
            }
        }
        BOOST_TEST_MESSAGE("Using CUDA device #" << cuda::driver::context::device());
    }

    ~set_cuda_device()
    {
        // Detach CUDA runtime from CUDA device context
        // This explicit clean-up is needed with CUDA < 3.0.
        cuda::thread::exit();
    }

private:
    boost::shared_ptr<cuda::driver::context> ctx_;
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
