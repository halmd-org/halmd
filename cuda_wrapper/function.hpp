/* cuda_wrapper/function.hpp
 *
 * Copyright (C) 2007  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef CUDA_FUNCTION_HPP
#define CUDA_FUNCTION_HPP

#include <cuda_wrapper/function_base.hpp>

namespace cuda
{

template <typename T>
class function;


template <typename T0>
class function<void (T0)> : public function_base
{
protected:
    typedef void T (T0);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


template <typename T0, typename T1>
class function<void (T0, T1)> : public function_base
{
protected:
    typedef void T (T0, T1);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


template <typename T0, typename T1, typename T2>
class function<void (T0, T1, T2)> : public function_base
{
protected:
    typedef void T (T0, T1, T2);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


template <typename T0, typename T1, typename T2, typename T3>
class function<void (T0, T1, T2, T3)> : public function_base
{
protected:
    typedef void T (T0, T1, T2, T3);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2, T3 x3)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	setup_argument(x3, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


template <typename T0, typename T1, typename T2, typename T3, typename T4>
class function<void (T0, T1, T2, T3, T4)> : public function_base
{
protected:
    typedef void T (T0, T1, T2, T3, T4);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2, T3 x3, T4 x4)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	setup_argument(x3, &offset);
	setup_argument(x4, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};

}

#endif /* ! CUDA_FUNCTION_HPP */
