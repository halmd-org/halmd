.. _download:

Downloading the source code
***************************

Cloning the repository
======================

HALMD is maintained in a public `Git <http://git-scm.com/>`_ repository ::

  git clone --recursive git://git.colberg.org/halmd.git

In case you are behind a firewall that blocks the git protocol port, use ::

  git clone --recursive http://git.colberg.org/halmd.git

This will create a directory ``halmd``, which holds a hidden copy of the
repository and a working copy of the source files.


Selecting a release version
===========================

By default, the above command yields the tip of the development branch.
A specific release version is checked out by ::

  git checkout TAG

where ``TAG`` is a valid release tag as listed from ::

  git tag -n3


Git tutorials
=============

If you are new to Git or version control in general, the `Git tutorial
<http://www.kernel.org/pub/software/scm/git/docs/gittutorial.html>`_
will get you started.

Former Subversion users may also read the `Git SVN Crash Course
<http://git.or.cz/course/svn.html>`_.

For in-depth documentation, see the `Git User's Manual
<http://www.kernel.org/pub/software/scm/git/docs/user-manual.html>`_.

