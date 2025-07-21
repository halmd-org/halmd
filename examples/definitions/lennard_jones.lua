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

-- define pair interaction with the (possibly truncated) Lennard-Jones potential,
-- use default parameters ε=1 and σ=1 for a single species
--
-- 'args' should provide entries for 'box', 'particle' and optionally 'cutoff' and 'smoothing'
--
-- returns the created force and potential modules
--

-- the following code is part of the documentation in doc/recipes/create_slit_pore.rst.in,
-- line number must not be changed
function M.create_pair_force(args)
    local box = utility.assert_kwarg(args, "box")
    local particle = utility.assert_kwarg(args, "particle")
    local cutoff = args.cutoff
    local smoothing = args.smoothing or 0.005

    -- define Lennard-Jones pair potential,
    local potential = mdsim.potentials.pair.lennard_jones()
    -- optionally, apply interaction cutoff
    if cutoff then
        -- use smooth truncation
        if smoothing > 0 then
            potential = potential:truncate({"smooth_r4", cutoff = cutoff, h = smoothing})
        else
            potential = potential:truncate({cutoff = cutoff})
        end
    end

    -- register computation of pair forces
    local force = mdsim.forces.pair({
        box = box, particle = particle, potential = potential
    })

    return force, potential
end

-- end of usage in doc/recipes/create_slit_pore.rst.in

return M
