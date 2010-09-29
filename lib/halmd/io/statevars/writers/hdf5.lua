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
local io = {
    statevars = {
        writers = {
            hdf5 = {
                [2] = halmd_wrapper.io.statevars.writers.hdf5_2_
              , [3] = halmd_wrapper.io.statevars.writers.hdf5_3_
            }
        }
    }
}
local args = require("halmd.options")
local assert = assert

module("halmd.io.statevars.writers.hdf5", halmd.modules.register)

--
-- construct HDF5 statevars writer module
--
function new()
    -- command line options
    local output = assert(args.output)
    local dimension = assert(args.dimension)

    -- parameters
    local file_name = output .. ".msv"

    return io.statevars.writers.hdf5[dimension](file_name)
end
