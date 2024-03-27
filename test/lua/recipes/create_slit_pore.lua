local mdsim = halmd.mdsim
local numeric = halmd.numeric

-- the following code is part of the documentation in doc/recipes/create_site_pore.rst.in,
-- line numbers must not be changed

--
-- define the particle-particle interaction using the truncated Lennard-Jones potential
--
function create_pair_force(box, particle, args)
    -- define pair potential
    local potential = mdsim.potentials.pair.lennard_jones()
    -- truncate potential
    if args.cutoff then
        if args.smoothing > 0 then
            potential = potential:truncate({"smooth_r4", cutoff = args.cutoff, h = args.smoothing})
        else
            potential = potential:truncate({cutoff = args.cutoff})
        end
    end

    return mdsim.forces.pair({box = box, particle = particle, potential = potential})
end

--
-- define the particle-wall interaction using the planar Lennard-Jones wall potential
--
function create_wall_force(box, particle, args)
    local wall_potential = mdsim.potentials.external.planar_wall({
        -- place two walls perpendicular to the x-axis
        surface_normal = {{-1, 0, 0}, {1, 0, 0}}
        -- add a space of .5σ between the nominal pore space and each of the walls
      , offset = {(args.pore_width + 1) / 2, (args.pore_width + 1) / 2}
        -- wall interaction parameters
      , sigma = 1, epsilon = args.epsilon, wetting = args.wetting
      , cutoff = math.pow(2, 1 / 6), smoothing = args.smoothing
      , species = 1
    })

    return mdsim.forces.external({box = box, particle = particle, potential = wall_potential})
end

function main()
local args = {
    box_length = 50, pore_width = 20
  , epsilon = 1, wetting = 1.5, smoothing = 0.005
  , density = 0.75, temperature  = 0.9
}

local length = {args.pore_width + 10, args.box_length, args.box_length}   -- add extra space "behind" the walls
local slab = {args.pore_width / length[1], 1, 1}
-- compute particle number from pore volume and mean density
local nparticle = math.floor(args.density * numeric.prod(slab) * numeric.prod(length))

-- create system state
local box = mdsim.box({length = length})
local particle = mdsim.particle({dimension = #length, particles = nparticle})

-- set initial particle positions
mdsim.positions.lattice({box = box, particle = particle, slab = slab}):set()

-- set initial particle velocities
mdsim.velocities.boltzmann({particle = particle, temperature = args.temperature}):set()

-- set up fluid-fluid forces
local ff_forces = create_pair_force(box, particle, args)

-- set up fluid-wall forces
local wall_forces = create_wall_force(box, particle, args)
-- end of usage in doc/usage/wetting.rst.in
end
