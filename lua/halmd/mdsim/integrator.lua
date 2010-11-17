--
-- Copyright © 2010  Peter Colberg and Felix Höfling
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
local integrator_wrapper = {
    [2] = halmd_wrapper.mdsim.integrator_2_
  , [3] = halmd_wrapper.mdsim.integrator_3_
  , nvt = {
      [2] = halmd_wrapper.mdsim.integrators.nvt_2_
    , [3] = halmd_wrapper.mdsim.integrators.nvt_3_
  }
}
local integrators = {
    verlet = require("halmd.mdsim.integrators.verlet")
  , verlet_nvt_andersen = require("halmd.mdsim.integrators.verlet_nvt_andersen")
}
local po = halmd_wrapper.po
local assert = assert

module("halmd.mdsim.integrator", halmd.modules.register)

--
-- construct integrator module
--
function new(args)
    local integrator = args.integrator or "verlet" --default value
    return integrators[integrator]()
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc, globals)
    globals:add("integrator", po.string(), "specify integration module")
end
