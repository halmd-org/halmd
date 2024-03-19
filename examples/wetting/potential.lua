#!/usr/bin/env halmd
--
-- Copyright © 2010-2023 Felix Höfling
-- Copyright © 2010-2012 Peter Colberg
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
-- Custom module for the definition of particle-particle and particle-wall interactions
--
local mdsim = require("halmd.mdsim")

local forces = {}

function forces.create_pair_force(box, particle, args)
-- define the particle-particle interaction using the truncated Lennard-Jones potential

    local lj_potential = mdsim.potentials.pair.lennard_jones()

    if args.cutoff > 0 then
        if args.smoothing > 0 then
                lj_potential = lj_potential:truncate({"smooth_r4", cutoff = args.cutoff, h = args.smoothing})
        else
                lj_potential = lj_potential:truncate({cutoff = args.cutoff})
        end
    end

    return mdsim.forces.pair({ box = box, particle = particle, potential = lj_potential})
end


function forces.create_wall_force(box, particle, args)
-- define the particle-wall interaction using the planar wall potential

    local wall_potential = mdsim.potentials.external.planar_wall({
        -- add space between particles and walls
        offset = {(args.pore_width + 1) / 2, (args.pore_width + 1) / 2} ,
        surface_normal = {{-1, 0, 0}, {1, 0,0 }},
        epsilon = args.epsilon,
        sigma = 1,
        wetting = args.wetting,
        cutoff = math.pow(2, 1/6),
        smoothing = args.smoothing,
        species = 1
    })

    return mdsim.forces.external({box = box, particle = particle, potential = wall_potential})
end

return forces