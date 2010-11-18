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

require("halmd.mdsim.integrators.verlet")
require("halmd.mdsim.integrators.verlet_nvt_andersen")

-- grab environment
local integrator_wrapper = {
    [2] = halmd_wrapper.mdsim.integrator_2_
  , [3] = halmd_wrapper.mdsim.integrator_3_
  , nvt = {
      [2] = halmd_wrapper.mdsim.integrators.nvt_2_
    , [3] = halmd_wrapper.mdsim.integrators.nvt_3_
  }
}
local integrators = halmd.mdsim.integrators
local po = halmd_wrapper.po
local assert = assert
local pairs = pairs

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

    -- integrator module choices with descriptions
    local choices = {}
    for integrator, module in pairs(integrators) do
        if module.name then
            choices[integrator] = module.name()
        end
    end

    globals:add("integrator", po.string():choices(choices), "specify integration module")
end
