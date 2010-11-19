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
local hdf5_writer = {
    host = {
        [2] = halmd_wrapper.io.trajectory.writers.hdf5_2_double_
      , [3] = halmd_wrapper.io.trajectory.writers.hdf5_3_double_
    }
  , gpu = {
        [2] = halmd_wrapper.io.trajectory.writers.hdf5_2_float_
      , [3] = halmd_wrapper.io.trajectory.writers.hdf5_3_float_
    }
}
local observables = {
    samples = {
        trajectory = require("halmd.observables.samples.trajectory")
    }
}
local mdsim = {
  core = require("halmd.mdsim.core")
}
local device = require("halmd.device")
local parameter = require("halmd.parameter")
local assert = assert

module("halmd.io.trajectory.writers.hdf5", halmd.modules.register)

--
-- construct HDF5 trajectory writer module
--
function new(args)
    -- dependency injection
    local core = mdsim.core()
    local dimension = assert(core.dimension)
    local sample = assert(observables.samples.trajectory())

    -- command line options
    local output = assert(args.output)

    -- parameters
    local file_name = output .. ".trj"

    local writer
    if not device() then
        writer = hdf5_writer.host[dimension](sample, file_name)
    else
        writer = hdf5_writer.gpu[dimension](sample, file_name)
    end

    parameter.register_writer(writer)

    return writer
end
