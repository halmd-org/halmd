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
-- definition of the modified Morse potential for fcc metals of Cu
--
-- Potential parameters were taken from:
-- MacDonald and MacDonald, Phys. Rev. B 24, 1715 (1981)
-- https://doi.org/10.1103/PhysRevB.24.1715
--

local mdsim = require("halmd.mdsim")
local utility = require("halmd.utility")
local constants = require("halmd.utility.constants")

local M = {}

-- define unit system: Ångstrøm, electron volt, picoseconds
-- conversion factors form simulation units to SI units
local eV = constants.eV
local Da = constants.Da
M.units = {
    length = 1e-10, time = 1e-12, energy = eV
  , mass = eV * math.pow(1e-12 / 1e-10, 2)      -- mass = energy × (time / length)^2
}

-- potential parameters in units of Ångstrøm and electron volt
-- atomic masses are taken from https://physics.nist.gov
M.parameters = {
    Cu = {"copper"
       , r_min = 2.5471
       , sigma = 1 / 1.1857
       , epsilon = 0.9403e-19 / eV
       , B = math.sqrt(2.265)
       , mass = 63.546 * Da / M.units.mass
    }
  , Ag = {"silver"
       , r_min = 2.8765
       , sigma = 1 / 1.1255
       , epsilon = 0.7874e-19 / eV
       , B = math.sqrt(2.3)
       , mass = 196.967 * Da / M.units.mass
    }
  , Ca = {"calcium"
       , r_min = 3.9264
       , sigma = 1 / 0.8380
       , epsilon = 0.5535e-19 / eV
       , B = math.sqrt(1.0)
       , mass = 40.078 * Da / M.units.mass
    }
  , Sr = {"strontium"
       , r_min = 4.2804
       , sigma = 1 / 0.7867
       , epsilon = 0.5442e-19 / eV
       , B = math.sqrt(1.0)
       , mass = 87.62 * Da / M.units.mass
    }
  , Al = {"aluminium"
       , r_min = 2.8485
       , sigma = 1 / 1.1611
       , epsilon = 0.6369e-19 / eV
       , B = math.sqrt(2.5)
       , mass = 26.982 * Da / M.units.mass
    }
  , Pb = {"lead"
       , r_min = 3.4779
       , sigma = 1 / 0.7776
       , epsilon = 0.5500e-19 / eV
       , B = math.sqrt(1.5)
       , mass = 207.2 * Da / M.units.mass
    }
  , Ni = {"nickel"
       , r_min = 2.4849
       , sigma = 1 / 1.3909
       , epsilon = 0.9843e-19 / eV
       , B = math.sqrt(2.4)
       , mass = 58.6934 * Da / M.units.mass
    }
}

--
-- set up pair interaction based on a modified Morse potential
--
-- 'args' must provide entries for 'box', 'particle', and 'substance',
-- optionally also 'cutoff' and 'smoothing' (note that the potential has been
-- developed without any cutoff). Additional keywords in 'args' are forwarded
-- to the constructor of the mdsim.forces.pair module.
--
-- possible values for 'substance' are the strings
-- "Cu", "Ag", "Ca", "Sr", "Al", "Pb", and "Ni".
--
-- returns the created force and potential modules
--
function M.create_pair_force(args)
    utility.assert_kwarg(args, "box")
    utility.assert_kwarg(args, "particle")

    local substance = utility.assert_kwarg(args, "substance")
    local param = M.parameters[substance]
    if not param then
        error(("unsupported chemical substance: %s"):format(substance), 3)
    end

    -- define Morse pair potential with parameters given above
    local potential = mdsim.potentials.pair.morse(param)

    -- optionally, apply interaction cutoff
    local cutoff = args.cutoff
    local smoothing = args.smoothing or 0.005
    if cutoff and cutoff > 0 then
        -- use smooth truncation
        if smoothing > 0 then
            potential = potential:truncate({"smooth_r4", cutoff = cutoff, h = smoothing})
        else
            potential = potential:truncate({cutoff = cutoff})
        end
    end

    -- register computation of pair forces,
    -- forward additional custom arguments to force module
    local force_args = args
    force_args.potential = potential
    local force = mdsim.forces.pair(force_args)

    return force, potential
end

return M
