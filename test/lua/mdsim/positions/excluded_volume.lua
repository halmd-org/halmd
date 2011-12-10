--
-- Copyright © 2011  Peter Colberg
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

require("halmd")

return function()
    local particle = halmd.mdsim.particle{
        particles = {5000}
      , masses = {1}
      , label = {"A"}
      , dimension = 3
    }
    halmd.mdsim.box{
        particles = {particle}
    }
    local sample = halmd.observables.samples.phase_space{
        memory = "host"
      , particle = particle
    }
    halmd.mdsim.positions.excluded_volume{
        samples = {sample}
      , diameters = {{1}}
    }
end
