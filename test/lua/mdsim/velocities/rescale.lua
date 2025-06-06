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

local function compute_and_assert_energy(particle, args)
    local mass = args.mass
    local stiffness = args.stiffness
    local dimension = #particle.data.velocity[1]
    local nparticle = #particle.data.velocity
    local expected = args.target_energy
    local tolerance = 1e-6

    local velocity = particle.data["velocity"]
    local position = particle.data["position"]

    local failed = 0
    for i = 1, nparticle do
        local v = velocity[i]
        local r = position[i]

        local v2, r2 = 0, 0
        for d = 1, dimension do
            v2 = v2 + v[d]^2
            r2 = r2 + r[d]^2
        end

        -- Check for invalid values
        if v2 ~= v2 or r2 ~= r2 then
            log.error(("NaN detected in particle %d: |v|² = %g, |r|² = %g"):format(i, v2, r2))
            failed = failed + 1
        else
            local kinetic = mass * v2 / 2
            local potential = stiffness * r2 / 2
            local energy = kinetic + potential

            assert(math.abs(energy - expected) <= tolerance * expected,
                ("particle #%d: difference (%g) between total energy (%g) and target value (%g) exceeds tolerance (%g)")
                   :format(i, math.abs(energy - expected), energy, expected, tolerance))
        end
    end
end

-- Shared test runner
local function run_rescale_test(args)
    local dimension = args.dimension
    local np = args.particles
    local L = args.box_length

    local length = {}
    for i = 1, dimension do length[i] = L end
    local box = mdsim.box({length = length})

    local particle = mdsim.particle({dimension = dimension, particles = np, species = 1})
    mdsim.positions.lattice({box = box, particle = particle}):set()

    local group = mdsim.particle_groups.all({particle = particle})

    if not args.zero_velocities then
        mdsim.velocities.boltzmann({
            particle = particle, group = group
          , temperature = args.temperature or 1
        }):set()
    end

    if args.stiffness > 0 then
        local stiffness = {}
        local offset = {}
        for i = 1, dimension do
            stiffness[i] = args.stiffness
            offset[i] = 0
        end
        local potential = mdsim.potentials.external.harmonic({
            stiffness = stiffness
          , offset = offset
        })
        mdsim.forces.external({box = box, particle = particle, potential = potential})
    end

    local thermo = observables.thermodynamics({group = group, box = box})

    local target_energy = args.target_energy
    local en_kin = thermo:kinetic_energy()
    local en_pot = thermo:potential_energy()

--    if target_energy <= en_pot then
--        error(("Invalid target energy: %.6f is less than potential energy %.6f"):format(target_energy, en_pot))
--    end

    log.info(("total energy:     %.6f"):format(en_kin + en_pot))
    log.info(("kinetic energy:   %.6f"):format(en_kin))
    log.info(("potential energy: %.6f"):format(en_pot))

    mdsim.velocities.rescale({
        particle = particle
      , target_energy = target_energy
    }):set()

    log.info(("total energy after rescaling:     %.6f"):format(thermo:internal_energy()))
    log.info(("kinetic energy after rescaling:   %.6f"):format(thermo:kinetic_energy()))
    log.info(("potential energy after rescaling: %.6f"):format(thermo:potential_energy()))

    compute_and_assert_energy(particle, args)
end

-- Test cases
test = {}

test["energy_target"] = function(args)
    args.target_energy = 50
    args.zero_velocities = false
    args.stiffness = 1
    run_rescale_test(args)
end

test["no_potential_energy"] = function(args)
    args.target_energy = 1
    args.zero_velocities = false
    args.stiffness = 0
    run_rescale_test(args)
end

test["zero_initial_velocity"] = function(args)
    args.target_energy = 50
    args.zero_velocities = true
    args.stiffness = 1
    run_rescale_test(args)
end

test["adding_mass"] = function(args)
    args.target_energy = 50
    args.zero_velocities = false
    args.stiffness = 1
    args.mass = 1
    run_rescale_test(args)
end


-- Argument parser
function define_args(parser)
    parser:add_argument("run_test", {
        type = "string",
        help = "select test: energy_target, no_potential_energy, zero_initial_velocity, adding_mass"
    })
    parser:add_argument("target_energy", {
        type = "number", default = 10,
        help = "target total energy to rescale velocities"
    })
    parser:add_argument("zero_velocities", {
        type = "boolean", default = false,
        help = "start with zero initial velocities"
    })
    parser:add_argument("stiffness", {
        type = "number", default = 1,
        help = "harmonic stiffness for potential energy"
    })
    parser:add_argument("mass", {type = "number", default = 1, help = "mass of particles"})
    parser:add_argument("particles", {type = "number", default = 100, help = "number of particles"})
    parser:add_argument("dimension", {type = "number", default = 3, help = "system dimensionality"})
    parser:add_argument("box-length", {type = "number", default = 10, help = "box's length"})
    parser:add_argument("temperature", {type = "number", default = 1.0, help = "temperature (if Boltzmann velocities)"})
end

-- Entry point
function main(args)
    -- run selected test case or, by default, all tests
    local test_case = args.run_test
    local cases = test_case and { test_case } or {"energy_target", "zero_initial_velocity", "no_potential_energy"}

    for i,case in ipairs(cases) do
        log.message(("Running test case '%s' ..."):format(case))
        assert(test[case], ("test case '%s' is not registered"):format(case))

        -- call test function
        test[case](args)
        log.message(("Test case '%s' finished."):format(case))
        log.message("")
    end
end
