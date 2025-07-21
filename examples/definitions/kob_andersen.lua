--
-- Copyright © 2023 Felix Höfling
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
-- custom module for the definition of particle interactions
--

local mdsim = require("halmd.mdsim")
local utility = require("halmd.utility")

local M = {}

-- define pair interaction for the binary Kob-Andersen mixture, which is based
-- on the (smoothly) truncated Lennard-Jones potential.
--
-- 'args' should provide entries for 'box', 'particle' and optionally 'smoothing'
--
-- returns the created force and potential modules
--
function M.create_pair_force(args)
    local box = utility.assert_kwarg(args, "box")
    local particle = utility.assert_kwarg(args, "particle")
    local cutoff = 2.5
    local smoothing = args.smoothing or 0.005

    -- define pair potential
    local potential = mdsim.potentials.pair.lennard_jones({
        epsilon = {
            {1  , 1.5} -- AA, AB
          , {1.5, 0.5} -- BA, BB
        }
      , sigma = {
            {1  , 0.8 } -- AA, AB
          , {0.8, 0.88} -- BA, BB
        }
    })
    -- apply smooth interaction cutoff
    if smoothing > 0 then
        potential = potential:truncate({"smooth_r4", cutoff = cutoff, h = smoothing})
    else
        potential = potential:truncate({cutoff = cutoff})
    end

    -- register computation of pair forces
    local force = mdsim.forces.pair({
        box = box, particle = particle, potential = potential
    })

    return force, potential
end

return M
