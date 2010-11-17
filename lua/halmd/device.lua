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
local device_wrapper
if halmd_wrapper.utility.gpu then
    device_wrapper = halmd_wrapper.utility.gpu.device
end
local po = halmd_wrapper.po
local assert = assert

module("halmd.device", halmd.modules.register)

local device -- singleton instance

--
-- construct device module
--
function new(args)
    if device_wrapper and not args.disable_gpu then
        local devices = args.devices or {}
        local threads = args.threads or 128 -- default value

        if not device then
            device = device_wrapper(devices, threads)
        end
        return device
    end
end

--
-- assemble module options
--
-- @param desc po.options_description
--
function options(desc)
    if device_wrapper then
        desc:add("devices,D", po.int_array(), "CUDA device(s)")
        desc:add("threads,T", po.uint(), "number of CUDA threads per block")
        desc:add("disable-gpu", po.bool_switch(), "disable GPU acceleration")
    end
end
