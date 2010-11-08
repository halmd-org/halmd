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
local hdf5_writer_wrapper = halmd_wrapper.io.profile.writers.hdf5
local args = require("halmd.options")
local parameter = require("halmd.parameter")
local assert = assert

module("halmd.io.profile.writers.hdf5", halmd.modules.register)

--
-- construct HDF5 profile writer module
--
function new()
    -- command line options
    local output = assert(args.output)

    local writer = hdf5_writer_wrapper(output .. ".prf")

    parameter.register_writer(writer)

    return writer
end
