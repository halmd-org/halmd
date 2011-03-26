Short guide on Luabind
----------------------

The *Luabind* library is used to export C++ classes to the Lua scripting
interface.

.. highlight:: c++

The typical wrapper pattern is illustrated with code from the Verlet integrator.
Let us start with the declaration of the abstract base class
``halmd::mdsim::integrator``. It is qualified for Lua export by the static
method ``luaopen``. ::

    // from file halmd/mdsim/integrator.hpp

    #include <lua.hpp>

    template <int dimension>
    class integrator
    {
    public:
        static void luaopen(lua_State* L);

        integrator() {}
        virtual ~integrator() {}
        virtual void integrate() = 0;
        virtual void finalize() = 0;
        virtual double timestep() const = 0;
        virtual void timestep(double timestep) = 0;
    };

The static method ``luaopen`` defines all entities that should be visible to Lua
by means of *Luabind*. It is called by the library constructor, i.e., before
programme execution starts. ::

    // from file halmd/mdsim/integrator.cpp

    #include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

    template <int dimension>
    void integrator<dimension>::luaopen(lua_State* L)
    {
        using namespace luabind;
        // dimension-dependent class name: integrator_2_, integrator_3_
        static std::string class_name("integrator_" + boost::lexical_cast<std::string>(dimension) + "_");
        // register a new Lua module
        module(L)
        [
            namespace_("halmd_wrapper")
            [
                namespace_("mdsim")
                [
                    class_<integrator, boost::shared_ptr<integrator> >(class_name.c_str())
                        // no constructor, this is an abstract base class
                        .property("timestep", (double (integrator::*)() const) &integrator::timestep)
                        .def("integrate", &integrator::integrate)
                        .def("finalize", &integrator::finalize)
                ]
            ]
        ];
    }

    HALMD_INIT( register_luaopen )
    {
        lua_wrapper::register_(0) //< distance of derived to base class
        [
            &integrator<3>::luaopen
        ]
        [
            &integrator<2>::luaopen
        ];
    }

**FIXME** explain the template arguments of ``class_``

**FIXME** explain the differences between property, def, def_readonly, add
comments in code sample

**FIXME** explain ``lua_wrapper::register_``

The actual class for the Verlet module derives from its abstract interface
class. Again, it has a static method ``luaopen``. Its constructor describes the
dependencies from other modules (particle, box) and specific parameters
(timestep). Further, it is good practice to define a type ``_Base`` pointing at
the base class.

The class contains another static method ``module_name``, which
is used by **FIXME** some Lua script to select the integrator via a global
command line option. ::

    // from file halmd/mdsim/host/integrators/verlet.hpp

    template <int dimension, typename float_type>
    class verlet
      : public mdsim::integrator<dimension>
    {
    public:
        typedef mdsim::integrator<dimension> _Base;
        typedef host::particle<dimension, float_type> particle_type;
        typedef mdsim::box<dimension> box_type;

        static char const* module_name() { return "verlet"; }

        boost::shared_ptr<particle_type> particle;
        boost::shared_ptr<box_type> box;

        static void luaopen(lua_State* L);

        verlet(
            boost::shared_ptr<particle_type> particle
          , boost::shared_ptr<box_type> box
          , double timestep
        );
        virtual void integrate();
        virtual void finalize();
        virtual void timestep(double timestep);
        virtual double timestep() const;
    };

Export to Lua is similar as for the base class. The main difference is that a
constructor is defined using ``def`` and that a wrapper is needed for the static
method ``module_name``. ::

    // from file halmd/mdsim/host/integrators/verlet.hpp

    template <int dimension, typename float_type>
    static char const* module_name_wrapper(verlet<dimension, float_type> const&)
    {
        return verlet<dimension, float_type>::module_name();
    }

    template <int dimension, typename float_type>
    void verlet<dimension, float_type>::luaopen(lua_State* L)
    {
        using namespace luabind;
        // dimension-dependent class name: verlet_2_, verlet_3_
        static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
        // register a new Lua module
        module(L)
        [
            namespace_("halmd_wrapper")
            [
                namespace_("mdsim")
                [
                    namespace_("host")
                    [
                        namespace_("integrators")
                        [
                            class_<verlet, shared_ptr<_Base>, bases<_Base> >(class_name.c_str())
                                .def(constructor<
                                    shared_ptr<particle_type>
                                  , shared_ptr<box_type>
                                  , double>()
                                )
                                .property("module_name", &module_name_wrapper<dimension, float_type>)
                        ]
                    ]
                ]
            ]
        ];
    }

    HALMD_INIT( register_luaopen )
    {
        lua_wrapper::register_(1) //< distance of derived to base class
        [
            &verlet<3, double>::luaopen
        ]
        [
            &verlet<2, double>::luaopen
        ];
    }

**FIXME** explain the three template arguments of ``class_``

**FIXME** explain why we need a wrapper for ``module_name``. Is it due to a
deficiency of luabind?

.. highlight:: lua

**FIXME** add some Lua code that exemplifies the usage of the exported module
::

    require("halmd.mdsim.integrator")
    integrator = assert(halmd_wrapper.mdsim.host.integrators.verlet_2_)
    print(integrator.module_name())
    instance = integrator(particle, box, 0.001)
    instance:integrate()
    instance:finalize()

**FIXME** when do you use a ``.`` and when a ``:`` for member access? Like ``core:run()`` but ``integrator.module_name()``?

