#
# Copyright © 2011  Peter Colberg
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# default installation prefix
PREFIX = $(HOME)/opt

# support parallel builds
ifdef CONCURRENCY_LEVEL
    PARALLEL_BUILD_FLAGS = -j$(CONCURRENCY_LEVEL)
endif

##
## define commonly used commands
##

WGET = wget -c
TAR = tar
RM = rm -rf
CP = cp -r
TOUCH = touch
PATCH = patch

##
## define top-level targets
##

build: build-cmake build-cmake build-lua build-boost build-luabind build-hdf5

fetch: fetch-cmake fetch-cmake fetch-lua fetch-boost fetch-luabind fetch-hdf5

install: install-cmake install-cmake install-lua install-boost install-luabind install-hdf5

clean: clean-cmake clean-cmake clean-lua clean-boost clean-luabind clean-hdf5

distclean: distclean-cmake distclean-cmake distclean-lua distclean-boost distclean-luabind distclean-hdf5

env: env-cmake env-cmake env-lua env-boost env-luabind env-hdf5

##
## CMake with CMake-CUDA patch
##

CMAKE_VERSION = 2.8.4
CMAKE_TARBALL = cmake-$(CMAKE_VERSION).tar.gz
CMAKE_TARBALL_URL = http://www.cmake.org/files/v2.8/$(CMAKE_TARBALL)
CMAKE_CUDA_PATCH = cmake-cuda-2.8.4-1-g2ed3c7a.patch
CMAKE_CUDA_PATCH_URL = http://sourceforge.net/projects/halmd/files/patches/$(CMAKE_CUDA_PATCH)
CMAKE_LIB64_PATCH = cmake_find_library_lib64_in_cxx_project.patch
CMAKE_LIB64_PATCH_URL = http://sourceforge.net/projects/halmd/files/patches/$(CMAKE_LIB64_PATCH)
CMAKE_BUILD_DIR = cmake-$(CMAKE_VERSION)
CMAKE_INSTALL_DIR = cmake-$(CMAKE_VERSION)

.fetch-cmake:
	$(WGET) $(CMAKE_TARBALL_URL)
	$(WGET) $(CMAKE_LIB64_PATCH_URL)
	$(WGET) $(CMAKE_CUDA_PATCH_URL)
	@$(TOUCH) $@

fetch-cmake: .fetch-cmake

.extract-cmake: .fetch-cmake
	$(RM) $(CMAKE_BUILD_DIR)
	$(TAR) -xzf $(CMAKE_TARBALL)
	cd $(CMAKE_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(CMAKE_CUDA_PATCH)
	cd $(CMAKE_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(CMAKE_LIB64_PATCH)
	@$(TOUCH) $@

extract-cmake: .extract-cmake

.configure-cmake: .extract-cmake
	cd $(CMAKE_BUILD_DIR) && ./configure --prefix=$(PREFIX)/$(CMAKE_INSTALL_DIR)
	@$(TOUCH) $@

configure-cmake: .configure-cmake

.build-cmake: .configure-cmake
	cd $(CMAKE_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-cmake: .build-cmake

install-cmake: .build-cmake
	cd $(CMAKE_BUILD_DIR) && make install

clean-cmake:
	@$(RM) .build-cmake
	@$(RM) .configure-cmake
	@$(RM) .extract-cmake
	$(RM) $(CMAKE_BUILD_DIR)

distclean-cmake: clean-cmake
	@$(RM) .fetch-cmake
	$(RM) $(CMAKE_TARBALL)
	$(RM) $(CMAKE_LIB64_PATCH)
	$(RM) $(CMAKE_CUDA_PATCH)

env-cmake:
	@echo
	@echo '# add CMake $(CMAKE_VERSION) to environment'
	@echo 'export PATH="$(PREFIX)/$(CMAKE_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(PREFIX)/$(CMAKE_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'

##
## Lua
##

LUA_VERSION = 5.1.4
LUA_TARBALL = lua-$(LUA_VERSION).tar.gz
LUA_TARBALL_URL = http://www.lua.org/ftp/$(LUA_TARBALL)
LUA_PATCH = patch-lua-5.1.4-2
LUA_PATCH_URL = http://www.lua.org/ftp/$(LUA_PATCH)
LUA_BUILD_DIR = lua-$(LUA_VERSION)
LUA_INSTALL_DIR = lua-$(LUA_VERSION)
LUA_CFLAGS = -DLUA_USE_LINUX -fPIC -O2 -Wall

.fetch-lua:
	$(WGET) $(LUA_TARBALL_URL)
	$(WGET) $(LUA_PATCH_URL)
	@$(TOUCH) $@

fetch-lua: .fetch-lua

.extract-lua: .fetch-lua
	$(RM) $(LUA_BUILD_DIR)
	$(TAR) -xzf $(LUA_TARBALL)
	cd $(LUA_BUILD_DIR)/src && $(PATCH) -p0 < $(CURDIR)/$(LUA_PATCH)
	@$(TOUCH) $@

extract-lua: .extract-lua

.build-lua: .extract-lua
	cd $(LUA_BUILD_DIR) && make linux CFLAGS="$(LUA_CFLAGS)"
	@$(TOUCH) $@

build-lua: .build-lua

install-lua: .build-lua
	cd $(LUA_BUILD_DIR) && make install INSTALL_TOP=$(PREFIX)/$(LUA_INSTALL_DIR)

clean-lua:
	@$(RM) .build-lua
	@$(RM) .extract-lua
	$(RM) $(LUA_BUILD_DIR)

distclean-lua: clean-lua
	@$(RM) .fetch-lua
	$(RM) $(LUA_TARBALL)
	$(RM) $(LUA_PATCH)

env-lua:
	@echo
	@echo '# add Lua $(LUA_VERSION) to environment'
	@echo 'export PATH="$(PREFIX)/$(LUA_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(PREFIX)/$(LUA_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(PREFIX)/$(LUA_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Boost C++ libraries with Boost.Log
##

BOOST_VERSION = 1.47.0
BOOST_RELEASE = 1_47_0
BOOST_TARBALL = boost_$(BOOST_RELEASE).tar.gz
BOOST_TARBALL_URL = http://sourceforge.net/projects/boost/files/boost/$(BOOST_VERSION)/$(BOOST_TARBALL)
BOOST_LOG_TARBALL = boost-log.tar
BOOST_LOG_TARBALL_URL = http://boost-log.svn.sourceforge.net/viewvc/boost-log/trunk/boost-log/?view=tar
BOOST_LOG_DIR = boost-log
BOOST_BUILD_DIR = boost_$(BOOST_RELEASE)
BOOST_INSTALL_DIR = boost_$(BOOST_RELEASE)
BOOST_CXXFLAGS = -fPIC

.fetch-boost:
	$(WGET) $(BOOST_TARBALL_URL)
	$(WGET) -O $(BOOST_LOG_TARBALL) $(BOOST_LOG_TARBALL_URL)
	@$(TOUCH) $@

fetch-boost: .fetch-boost

.extract-boost: .fetch-boost
	$(RM) $(BOOST_BUILD_DIR) $(BOOST_LOG_DIR)
	$(TAR) -xzf $(BOOST_TARBALL)
	$(TAR) -xf $(BOOST_LOG_TARBALL)
	$(CP) $(BOOST_LOG_DIR)/boost/log $(BOOST_BUILD_DIR)/boost/
	$(CP) $(BOOST_LOG_DIR)/libs/log $(BOOST_BUILD_DIR)/libs/
	@$(TOUCH) $@

extract-boost: .extract-boost

.configure-boost: .extract-boost
	cd $(BOOST_BUILD_DIR) && ./bootstrap.sh
	@$(TOUCH) $@

configure-boost: .configure-boost

.build-boost: .configure-boost
	cd $(BOOST_BUILD_DIR) && ./bjam cxxflags="$(BOOST_CXXFLAGS)" $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-boost: .build-boost

install-boost: .build-boost
	cd $(BOOST_BUILD_DIR) && ./bjam cxxflags="$(BOOST_CXXFLAGS)" install --prefix=$(PREFIX)/$(BOOST_INSTALL_DIR)

clean-boost:
	@$(RM) .build-boost
	@$(RM) .configure-boost
	@$(RM) .extract-boost
	$(RM) $(BOOST_BUILD_DIR) $(BOOST_LOG_DIR)

distclean-boost: clean-boost
	@$(RM) .fetch-boost
	$(RM) $(BOOST_TARBALL)
	$(RM) $(BOOST_LOG_TARBALL)

env-boost:
	@echo
	@echo '# add Boost $(BOOST_VERSION) to environment'
	@echo 'export PATH="$(PREFIX)/$(BOOST_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(PREFIX)/$(BOOST_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export PYTHONPATH="$(PREFIX)/$(BOOST_INSTALL_DIR)/lib$${PYTHONPATH+:$$PYTHONPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(PREFIX)/$(BOOST_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Luabind library with Clang C++ compiler fix
##

LUABIND_VERSION = 0.9.1
LUABIND_TARBALL = luabind-$(LUABIND_VERSION).tar.gz
LUABIND_TARBALL_URL = http://sourceforge.net/projects/luabind/files/luabind/$(LUABIND_VERSION)/$(LUABIND_TARBALL)
LUABIND_PATCH = luabind_proper_forward_declarations_for_clang.patch
LUABIND_PATCH_URL = http://sourceforge.net/projects/halmd/files/patches/$(LUABIND_PATCH)
LUABIND_BUILD_DIR = luabind-$(LUABIND_VERSION)
LUABIND_BUILD_FLAGS = cxxflags=-fPIC link=static variant=release variant=debug
LUABIND_INSTALL_DIR = luabind-$(LUABIND_VERSION)

.fetch-luabind:
	$(WGET) $(LUABIND_TARBALL_URL)
	$(WGET) $(LUABIND_PATCH_URL)
	@$(TOUCH) $@

fetch-luabind: .fetch-luabind

.extract-luabind: .fetch-luabind .build-lua
	$(RM) $(LUABIND_BUILD_DIR)
	$(TAR) -xzf $(LUABIND_TARBALL)
	cd $(LUABIND_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(LUABIND_PATCH)
	mkdir $(LUABIND_BUILD_DIR)/lua
	ln -s $(CURDIR)/$(LUA_BUILD_DIR)/src $(LUABIND_BUILD_DIR)/lua/include
	ln -s $(CURDIR)/$(LUA_BUILD_DIR)/src $(LUABIND_BUILD_DIR)/lua/lib
	@$(TOUCH) $@

extract-luabind: .extract-luabind

.build-luabind: .extract-luabind .configure-boost
	cd $(LUABIND_BUILD_DIR) && BOOST_ROOT=$(CURDIR)/$(BOOST_BUILD_DIR) LUA_PATH=$(CURDIR)/$(LUABIND_BUILD_DIR)/lua $(CURDIR)/$(BOOST_BUILD_DIR)/bjam $(LUABIND_BUILD_FLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-luabind: .build-luabind

install-luabind: .build-luabind
	cd $(LUABIND_BUILD_DIR) && BOOST_ROOT=$(CURDIR)/$(BOOST_BUILD_DIR) LUA_PATH=$(CURDIR)/$(LUABIND_BUILD_DIR)/lua $(CURDIR)/$(BOOST_BUILD_DIR)/bjam $(LUABIND_BUILD_FLAGS) install --prefix=$(PREFIX)/$(LUABIND_INSTALL_DIR)

clean-luabind:
	@$(RM) .build-luabind
	@$(RM) .extract-luabind
	$(RM) $(LUABIND_BUILD_DIR)

distclean-luabind: clean-luabind
	@$(RM) .fetch-luabind
	$(RM) $(LUABIND_TARBALL)
	$(RM) $(LUABIND_PATCH)

env-luabind:
	@echo
	@echo '# add Luabind $(LUABIND_VERSION) to environment'
	@echo 'export LD_LIBRARY_PATH="$(PREFIX)/$(LUABIND_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(PREFIX)/$(LUABIND_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## HDF5 C++ library
##

HDF5_VERSION = 1.8.7
HDF5_TARBALL = hdf5-$(HDF5_VERSION).tar.bz2
HDF5_TARBALL_URL = http://www.hdfgroup.org/ftp/HDF5/current/src/$(HDF5_TARBALL)
HDF5_BUILD_DIR = hdf5-$(HDF5_VERSION)
HDF5_INSTALL_DIR = hdf5-$(HDF5_VERSION)
HDF5_CONFIGURE_FLAGS = --enable-cxx --enable-static --disable-shared
HDF5_CFLAGS = -fPIC
HDF5_CXXFLAGS = -fPIC

.fetch-hdf5:
	$(WGET) $(HDF5_TARBALL_URL)
	@$(TOUCH) $@

fetch-hdf5: .fetch-hdf5

.extract-hdf5: .fetch-hdf5
	$(RM) $(HDF5_BUILD_DIR)
	$(TAR) -xjf $(HDF5_TARBALL)
	@$(TOUCH) $@

extract-hdf5: .extract-hdf5

.configure-hdf5: .extract-hdf5
	cd $(HDF5_BUILD_DIR) && CFLAGS="$(HDF5_CFLAGS)" CXXFLAGS="$(HDF5_CXXFLAGS)" ./configure $(HDF5_CONFIGURE_FLAGS) --prefix=$(PREFIX)/$(HDF5_INSTALL_DIR)
	@$(TOUCH) $@

configure-hdf5: .configure-hdf5

.build-hdf5: .configure-hdf5
	cd $(HDF5_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-hdf5: .build-hdf5

install-hdf5: .build-hdf5
	cd $(HDF5_BUILD_DIR) && make install

clean-hdf5:
	@$(RM) .build-hdf5
	@$(RM) .configure-hdf5
	@$(RM) .extract-hdf5
	$(RM) $(HDF5_BUILD_DIR)

distclean-hdf5: clean-hdf5
	@$(RM) .fetch-hdf5
	$(RM) $(HDF5_TARBALL)

env-hdf5:
	@echo
	@echo '# add HDF5 $(HDF5_VERSION) to environment'
	@echo 'export PATH="$(PREFIX)/$(HDF5_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(PREFIX)/$(HDF5_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(PREFIX)/$(HDF5_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Git version control
##

GIT_VERSION = 1.7.6
GIT_TARBALL = git-$(GIT_VERSION).tar.bz2
GIT_MANPAGES_TARBALL = git-manpages-$(GIT_VERSION).tar.bz2
GIT_TARBALL_URL = http://kernel.org/pub/software/scm/git/$(GIT_TARBALL)
GIT_MANPAGES_TARBALL_URL = http://kernel.org/pub/software/scm/git/$(GIT_MANPAGES_TARBALL)
GIT_BUILD_DIR = git-$(GIT_VERSION)
GIT_INSTALL_DIR = git-$(GIT_VERSION)

.fetch-git:
	$(WGET) $(GIT_TARBALL_URL)
	$(WGET) $(GIT_MANPAGES_TARBALL_URL)
	@$(TOUCH) $@

fetch-git: .fetch-git

.extract-git: .fetch-git
	$(RM) $(GIT_BUILD_DIR)
	$(TAR) -xjf $(GIT_TARBALL)
	@$(TOUCH) $@

extract-git: .extract-git

.configure-git: .extract-git
	cd $(GIT_BUILD_DIR) && ./configure --prefix=$(PREFIX)/$(GIT_INSTALL_DIR)
	@$(TOUCH) $@

configure-git: .configure-git

.build-git: .configure-git
	cd $(GIT_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-git: .build-git

install-git: .build-git
	cd $(GIT_BUILD_DIR) && make install
	install -d $(PREFIX)/$(GIT_INSTALL_DIR)/share/man
	cd $(PREFIX)/$(GIT_INSTALL_DIR)/share/man && $(TAR) -xjf $(CURDIR)/$(GIT_MANPAGES_TARBALL)

clean-git:
	@$(RM) .build-git
	@$(RM) .configure-git
	@$(RM) .extract-git
	$(RM) $(GIT_BUILD_DIR)

distclean-git: clean-git
	@$(RM) .fetch-git
	$(RM) $(GIT_TARBALL)
	$(RM) $(GIT_MANPAGES_TARBALL)

env-git:
	@echo
	@echo '# add Git $(GIT_VERSION) to environment'
	@echo 'export PATH="$(PREFIX)/$(GIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(PREFIX)/$(GIT_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'
