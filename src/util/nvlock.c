/*
 * nvlock - Exclusively lock an unused NVIDIA device and execute given program
 *
 * Copyright (C) 2008  Peter Colberg
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

#include <fcntl.h>
#include <glob.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char** argv)
{
    glob_t globbuf;
    int i, fd;

    if (argc < 2) {
	fprintf(stderr, "usage: nvlock <command>\n");
	return 1;
    }

    glob("/dev/nvidia[0-9]", GLOB_DOOFFS, NULL, &globbuf);

    if (!globbuf.gl_pathc) {
	fprintf(stderr, "nvlock: no NVIDIA devices found\n");
	return 1;
    }

    /* use devices in reverse order as first device is CUDA default */
    for (i = globbuf.gl_pathc - 1; i >= 0; --i) {
	if (-1 == (fd = open(globbuf.gl_pathv[i], O_RDWR))) {
            fprintf(stderr, "nvlock: could not open device: %s", globbuf.gl_pathv[i]);
	    return 1;
	}
	if (-1 == flock(fd, LOCK_EX | LOCK_NB)) {
	    close(fd);
	    continue;
	}

	/* set numeric CUDA device in environment variable */
        if (-1 == setenv("CUDA_DEVICE", globbuf.gl_pathv[i][11], 1)) {
	    fprintf(stderr, "nvlock: failed to set CUDA_DEVICE environment variable\n");
	    return 1;
	}
        /* replace process with given program */
        if (-1 == execvp(argv[1], &argv[1])) {
	    fprintf(stderr, "nvlock: failed to execute process: %s\n", argv[1]);
	    return 1;
	}
    }

    fprintf(stderr, "nvlock: no unused NVIDIA device found\n");
    return 1;
}
