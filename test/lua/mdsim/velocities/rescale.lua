--
-- Copyright © 2025 Giorgia Marcelli
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

local log = halmd.io.log
local mdsim = halmd.mdsim
local observables = halmd.observables

-- Shared setup function
local function run_rescale_test(args)
    local dimension = args.dimension
    local np = args.particles
    local L = args.box_length

    local box_length = {}
    for i = 1, dimension do box_length[i] = L end
    local box = mdsim.box({length = box_length})

    -- Particle system
    local particle = mdsim.particle({dimension = dimension, particles = np, species = 1})
    mdsim.positions.lattice({box = box, particle = particle}):set()

    local group = mdsim.particle_groups.all({particle = particle})

    -- Velocity assignment
    if not args.zero_velocities then
        mdsim.velocities.boltzmann({
            particle = particle,
            group = group,
            temperature = args.temperature or 1.0
        }):set()
    end

    -- Optional harmonic potential
    if args.stiffness > 0 then
        local potential = mdsim.potentials.external.harmonic({
            stiffness = {args.stiffness},
            offset = {0}
        })
        mdsim.forces.external({box = box, particle = particle, potential = potential})
    end

    -- Thermodynamics and rescaling
    local thermo = observables.thermodynamics({box = box, group = group})

    mdsim.velocities.rescale({
        particle = particle,
        energy = args.target_energy
    }):set()

    -- Logging
    log.message(("Target energy: %.3f"):format(args.target_energy))
    log.message(("Kinetic energy: %.6f"):format(thermo.kinetic_energy()))
    log.message(("Potential energy: %.6f"):format(thermo.potential_energy()))
    log.message(("Total energy: %.6f"):format(thermo.energy()))
end

-- Test cases
test = {}

-- 1. E=1, Boltzmann, stiffness=1
test["Etot=1"] = function(args)
    args.target_energy = 1
    args.zero_velocities = false
    args.stiffness = 1
    run_rescale_test(args)
end

-- 2. E=2, Boltzmann, stiffness=1
test["Etot=2"] = function(args)
    args.target_energy = 2
    args.zero_velocities = false
    args.stiffness = 1
    run_rescale_test(args)
end

-- 3. E=1, Boltzmann, stiffness=0
test["K=0"] = function(args)
    args.target_energy = 1
    args.zero_velocities = false
    args.stiffness = 0
    run_rescale_test(args)
end

-- 4. E=1, Zero velocities, stiffness=1
test["v=0"] = function(args)
    args.target_energy = 1
    args.zero_velocities = true
    args.stiffness = 1
    run_rescale_test(args)
end

-- Argument parsing
function define_args(parser)
    parser:add_argument("run_test", {type = "string", help = "Select case1, case2, case3, or case4"})
    parser:add_argument("particles", {type = "number", default = 10000, help = "Number of particles"})
    parser:add_argument("dimension", {type = "number", default = 3, help = "Simulation dimension"})
    parser:add_argument("box-length", {type = "number", default = 10, help = "Box side length"})
    parser:add_argument("temperature", {type = "number", default = 1.0, help = "Temperature for Boltzmann velocities"})
end

-- Run selected case
function main(args)
    assert(test[args.run_test], "Unknown test: " .. args.run_test)
    log.message(("Running test: %s"):format(args.run_test))
    test[args.run_test](args)
    log.message(("Finished test: %s"):format(args.run_test))
end
