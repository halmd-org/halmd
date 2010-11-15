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
local force_wrapper = {
    [2] = halmd_wrapper.mdsim.force_2_
  , [3] = halmd_wrapper.mdsim.force_3_
}
local forces = {
    pair_trunc = require("halmd.mdsim.forces.pair_trunc")
}
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.force", halmd.modules.register)

--
-- construct force module
--
function new(args)
    local force = args.force or "lennard_jones" -- default value
    return forces.pair_trunc{ force = force }
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    desc:add("force", po.string(), "specify force module")
end
