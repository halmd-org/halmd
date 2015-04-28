--
-- Copyright © 2013      Nicolas Höft
-- Copyright © 2010-2015 Felix Höfling
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
-- custom module for the rescaling of particle velocities after thermalisation
--
local M = {}

local log = require("halmd.io.log")
local utility = require("halmd.utility")
local assert = assert
local error = error

--
-- Shift velocities to zero mean and rescale them to the given internal energy
-- or temperature
--
-- :param table args: keyword arguments
-- :param args.msv: instance of :mod:`halmd.observables.thermodynamics`
-- :param number args.internal_energy: target value for internal energy
-- :param number args.temperature: target value for (instantaneous) temperature
--
-- :returns: applied rescaling factor
--
local function rescale_velocity(self, args)
    local msv   = utility.assert_kwarg(args, "msv")
    local group = assert(msv.group)
    local internal_energy = args.internal_energy
    local temperature = args.temperature
    if not internal_energy and not temperature then
        error("missing keyword argument: temperature or internal_energy", 2)
    end

    -- determine centre-of-mass velocity and subtract it
    local vcm  = msv:center_of_mass_velocity()
    local dimension = #vcm
    for d = 1, dimension do
        vcm[d] = -vcm[d]
    end
    group.particle:shift_velocity_group(group, vcm)

    -- compute kinetic energy assuming vcm=0
    local ekin = msv:kinetic_energy()

    local scale = 1
    if internal_energy then
        -- internal energy = potential energy + kinetic energy
        local epot = msv:potential_energy()
        if(internal_energy < epot) then
            error("velocity rescaling failed: internal energy too small")
        end
        scale = math.sqrt((internal_energy - epot) / ekin)
        log.info("target value of internal energy: %f", internal_energy)
    elseif temperature then
        -- rescale velocities to the prescribed temperature (kinectic energy)
        scale = math.sqrt(temperature * dimension / ekin)
        log.info("target temperature: %f", temperature)
    end

    log.info("rescale velocities of %s particles by a factor of %f", group.label, scale)
    group.particle:rescale_velocity_group(group, scale)
    return scale
end

-- override the __call metamethod so that we can use the module as a functor
setmetatable(M, {__call = rescale_velocity})

return M
