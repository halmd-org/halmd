#!/usr/bin/env halmd
--
-- Copyright © 2020 Felix Höfling
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as
-- published by the Free Software Foundation, either version 3 of
-- the License, or (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General
-- Public License along with this program.  If not, see
-- <http://www.gnu.org/licenses/>.
--

-- grab modules
local log = halmd.io.log

function string_interp()
    local input = "{name:s} is {val:7.2f}%"
    local output = input:interp({name = "concentration_NOx", val = 56.2795})

    log.message("formatting string: %s", input)
    log.message("interpolated string: %s", output)
    assert(output == "concentration_NOx is   56.28%")
end

function main()
    string_interp()
end
