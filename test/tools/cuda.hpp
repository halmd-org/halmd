/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_TEST_TOOLS_CUDA_HPP
#define HALMD_TEST_TOOLS_CUDA_HPP

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_monitor.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <test/tools/init.hpp>

/**
 * "global fixture" for Boost Unit Test Framework: select CUDA device
 *
 * create a CUDA context on a free device, respect device locks
 */
struct set_cuda_device
{
    /**
     * implicitly create CUDA context
     */
    set_cuda_device()
    {
        cuda::thread::synchronize();
        BOOST_TEST_MESSAGE( "Using CUDA device #" << cuda::device::get() );
    }

    /**
     * destroy CUDA context
     */
    ~set_cuda_device()
    {
        cuda::thread::exit();
    }
};

void cuda_error_translator(cuda::error const& e)
{
    BOOST_TEST_MESSAGE( "(CUDA error) " << e.what() );
    throw;
}

HALMD_TEST_INIT( register_cuda_error )
{
    using namespace boost::unit_test;
    unit_test_monitor.register_exception_translator<cuda::error>(&cuda_error_translator);
}

#endif /* ! HALMD_TEST_TOOLS_CUDA_HPP */
