#!/usr/bin/env halmd
--
-- Copyright © 2014 Nicolas Höft
-- Copyright © 2014 Felix Höfling
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.
--

local halmd = require("halmd")
halmd.io.log.open_console()

-- grab modules
local log = halmd.io.log
local numeric = halmd.utility.numeric

function test(dims)
    local ndim = 1
    for d = 1, #dims do
        ndim = ndim * dims[d]
    end
    log.info("total number of points: %d", ndim)

    for i = 1, ndim do
        local index = numeric.offset_to_multi_index(i, dims)
        local offset = numeric.multi_index_to_offset(index, dims)
        log.info("offset: %d, index: %d %d %d %d", offset, index[1], index[2], index[3], index[4])
        assert(i == offset)
    end
end

test({2, 1, 3, 3})
test({1})
