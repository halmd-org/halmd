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

--
-- HALMD Lua library
--

require("halmd.gpu.device")
require("halmd.gpu.random")
require("halmd.io.logging")
require("halmd.io.trajectory.reader")
require("halmd.mdsim.box")
require("halmd.mdsim.core")
require("halmd.mdsim.force")
require("halmd.mdsim.gpu.neighbour")
require("halmd.mdsim.host.forces.lj")
require("halmd.mdsim.host.forces.power_law")
require("halmd.mdsim.host.forces.smooth")
require("halmd.mdsim.host.velocities.boltzmann")
require("halmd.mdsim.host.neighbour")
require("halmd.mdsim.integrator")
require("halmd.mdsim.particle")
require("halmd.mdsim.position")
require("halmd.mdsim.velocity")
require("halmd.observables.thermodynamics")
require("halmd.random")
require("halmd.sampler")

-- grab environment
local log = log

module("halmd")

function run()
    log.info("********* Hello, World! *********")
end
