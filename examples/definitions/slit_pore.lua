#!/usr/bin/env halmd
--
-- Copyright © 2024 Felix Höfling
-- Copyright © 2024 Jung Nguyen
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as
-- published by the Free Software Foundation, either version 3 of
-- the License, or (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General
-- Public License along with this program.  If not, see
-- <http://www.gnu.org/licenses/>.
--

--
-- Custom module for the definition of particle–wall interactions
-- realising confinement in a slit pore
--
local mdsim = require("halmd.mdsim")
local utility = require("halmd.utility")

local M = {}

--
-- define the particle-wall interaction using the planar Lennard-Jones wall potential
--
-- 'args' should provide entries for 'box', 'particle', 'pore_width', and the interaction parameters
--
-- returns the created force and potential modules

-- the following code is part of the documentation in doc/recipes/create_slit_pore.rst.in,
-- line number must not be changed
function M.create_wall_force(args)
    local box = utility.assert_kwarg(args, "box")
    local particle = utility.assert_kwarg(args, "particle")
    local pore_width = utility.assert_kwarg(args, "pore_width")
    local sigma = args.sigma or 1
    local epsilon = args.epsilon or 1
    local wetting = args.wetting or 1
    local cutoff = args.cutoff or 2.5
    local smoothing = args.smoothing or 0.005

    local potential = mdsim.potentials.external.planar_wall({
        -- place two walls perpendicular to the x-axis
        surface_normal = {{-1, 0, 0}, {1, 0, 0}}
        -- add a space of .5σ between the nominal pore space and each of the walls
      , offset = {(pore_width + 1) / 2, (pore_width + 1) / 2}
        -- wall interaction parameters
      , sigma = sigma, epsilon = epsilon, wetting = wetting
      , cutoff = cutoff,  smoothing = smoothing
      , species = (type(sigma) == 'number') and 1 or nil
    })

    -- register computation of wall forces
    local force = mdsim.forces.external({
        box = box, particle = particle, potential = potential
    })

    return force, potential
end

-- end of usage in doc/recipes/create_slit_pore.rst.in

return M
