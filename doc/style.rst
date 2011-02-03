Coding conventions
------------------

#) New modules undergo rigorous peer reviewing before entering the branch of a release candidate.

#) Code design should observe the guidelines of H. Sutter and A. Alexandrescu in "C++ Coding Standards" (Addison Wesely, 2005).

#) TODO examplatory code fragments, formatting guidelines

    .. code-block:: c

        for (unsigned int i = 0; i < count; ++i) {    //< use pre-increment
            // loop body
        }

#) Naming rules

    * all identifiers are in lower case throughout, long names may be grouped by underscore, type identifiers end in ``_type``

    * private class attributes and methods end with an underscore

    * names of containers and collections are singular

#) Use exception-safe containers and pointers

    .. todo::

       RAII, STL containers, boost::shared_ptr
