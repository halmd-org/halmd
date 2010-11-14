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
local velocities = {
    boltzmann = require("halmd.mdsim.velocities.boltzmann")
}
local velocity_wrapper = {
    [2] = halmd_wrapper.mdsim.velocity_2_
  , [3] = halmd_wrapper.mdsim.velocity_3_
}
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.velocity", halmd.modules.register)

function options(desc)
    desc:add("velocity", po.string(), "initial particle velocities module")
end

--
-- construct velocity module
--
function new(args)
    local velocity = args.velocity or "boltzmann" -- default value
    return velocities[velocity]()
end
