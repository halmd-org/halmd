/*
 * nvlock - Exclusively lock an unused NVIDIA device and execute given program
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

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    if (argc < 2) {
	fprintf(stderr, "usage: nvlock <command>\n");
	return 1;
    }
    if (-1 == setenv("LD_PRELOAD", "libnvlock.so", 1)) {
	fprintf(stderr, "nvlock: failed to set LD_PRELOAD environment variable\n");
	return 1;
    }
    if (-1 == execvp(argv[1], &argv[1])) {
	fprintf(stderr, "nvlock: failed to execute process: %s\n", argv[1]);
	return 1;
    }
    return 0;
}
