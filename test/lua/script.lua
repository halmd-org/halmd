-- the following code is part of the documentation in doc/usage/script.rst.in
-- grab module namespaces
local log = halmd.io.log
local random = halmd.random

function main(args)
    -- some output to logger
    log.info("Write 'Hello World!' to " .. args.output .. ".log")

    -- seed the random number generator
    if args.random_seed then
      random.generator({seed = args.random_seed})
    end

    -- here: setup system and run simulation
end

function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.substitute_date_time,
        default = "project_%Y%m%d_%H%M%S", help = "prefix of output files"})

    parser:add_argument("random-seed", {type = "integer", help = "seed for random number generator"})

    return parser
end
-- end of usage in doc/recipes/create_mixture.rst.in
