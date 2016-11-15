function run()
-- the following code is part of the documentation in doc/recipes/create_mixture.rst.in

local mdsim   = halmd.mdsim
local numeric = halmd.numeric
local random  = halmd.random

local nparticle = {8000, 2000} -- particle numbers for each component
local length = {20, 20, 20}    -- cubic simulation box
local temperature = 1.5        -- temperature of Maxwell–Boltzmann distribution

-- create system state
local box = mdsim.box({length = length})
local particle = mdsim.particle({dimension = #length, particles = numeric.sum(nparticle), species = #nparticle})

-- assign particle species
local species = {}
for s = 1, #nparticle do
    for i = 1, nparticle[s] do
        table.insert(species, s - 1)    -- species values are 0-based 
    end
end
particle:set_species(species)

-- set particle positions sequentially on an fcc lattice
mdsim.positions.lattice({box = box, particle = particle}):set()

-- shuffle positions randomly
local r = particle:get_position()
r = random.generator({memory = "host"}):shuffle(r)
particle:set_position(r)

-- set initial particle velocities
mdsim.velocities.boltzmann({particle = particle, temperature = temperature}):set()

-- end of usage in doc/recipes/create_mixture.rst.in
end