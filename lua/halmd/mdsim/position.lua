--
-- Copyright © 2010  Peter Colberg
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

require("halmd.modules")

-- grab environment
local position_wrapper = {
    [2] = halmd_wrapper.mdsim.position_2_
  , [3] = halmd_wrapper.mdsim.position_3_
}
local positions = {
    lattice = require("halmd.mdsim.position.lattice")
}
local assert = assert

module("halmd.mdsim.position", halmd.modules.register)

--
-- construct position module
--
function new(args)
    local position = assert(args.position)
    return positions[position]()
end

options = position_wrapper[2].options

