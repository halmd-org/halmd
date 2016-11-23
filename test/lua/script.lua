-- the following code is part of the documentation in doc/usage/script.rst.in
-- grab module namespaces
local log = halmd.io.log

function main(args)
    -- some output to logger
    log.info("Write 'Hello World!' to " .. args.output .. ".log")

    -- here: setup system and run simulation
end

function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.substitute_date_time_action,
        default = "project_%Y%m%d_%H%M%S", help = "prefix of output files"})

    parser:add_argument("random-seed", {type = "integer", action = parser.random_seed_action,
        help = "seed for random number generator"})

    return parser
end
-- end of usage in doc/recipes/create_mixture.rst.in
