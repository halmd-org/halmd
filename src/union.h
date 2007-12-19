/* union.h
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

#ifndef __UNION_H__
#define __UNION_H__


#define __offsetof(var, member) offsetof(__typeof(var), member)


template <typename X, typename Y>
struct union2
{
  X x; Y y;
};


template <typename X, typename Y, typename Z>
struct union3
{
  X x; Y y; Z z;
};


template <typename X, typename Y, typename Z, typename W>
struct union4
{
  X x; Y y; Z z; W w;
};


#endif /* ! __UNION_H__ */
