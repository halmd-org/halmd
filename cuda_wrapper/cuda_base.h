/* cuda_base.h
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

#ifndef __CUDA_BASE_H__
#define __CUDA_BASE_H__

#include <cuda_runtime.h>


class cuda_base
{
public:
  cuda_base()
  {
  }

private:
  /* disable default copy constructor */
  cuda_base(const cuda_base&);
  /* disable default assignment operator */
  void operator=(const cuda_base&);
  /* disable default address-of operator */
  void *operator&();
};


#endif /* ! __CUDA_BASE_H__ */
