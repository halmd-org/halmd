#
# Copyright © 2011-2013 Peter Colberg
# Copyright © 2013      Nicolas Höft
# Copyright © 2013      Felix Höfling
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

BASE64    = base64 -di
CMAKE     = cmake
CP        = cp -r
GIT       = git
GUNZIP    = gzip -d
MKDIR     = mkdir -p
PATCH     = patch
RM        = rm -rf
SED       = sed
SHA256SUM = sha256sum --check
TAR       = tar
TOUCH     = touch
WGET      = wget

##
## define top-level targets
##

build: build-cmake build-lua build-boost build-hdf5

fetch: fetch-cmake fetch-lua fetch-boost fetch-hdf5

install: install-cmake install-lua install-boost install-hdf5

clean: clean-cmake clean-lua clean-boost clean-hdf5

distclean: distclean-cmake distclean-lua distclean-boost distclean-hdf5

env: env-cmake env-lua env-boost env-hdf5

##
## CMake with CMake-CUDA patch
##

CMAKE_CUDA_VERSION = 2.8.11.2
CMAKE_CUDA_TARBALL = $(CMAKE_CUDA_VERSION).tar.gz
CMAKE_CUDA_TARBALL_SHA256 = 1a498c187bd038ceeb7f0409fb59d849b74871f4c0f06301999d75872175a5b9
CMAKE_CUDA_TARBALL_URL = http://git.halmd.org/cmake-cuda/archive/$(CMAKE_CUDA_TARBALL)
CMAKE_SOURCE_DIR = cmake-cuda-$(CMAKE_CUDA_VERSION)
CMAKE_BUILD_DIR = $(CMAKE_SOURCE_DIR)/build
CMAKE_INSTALL_DIR = $(PREFIX)/cmake-$(CMAKE_CUDA_VERSION)

.fetch-cmake-$(CMAKE_CUDA_VERSION):
	@$(RM) $(CMAKE_CUDA_TARBALL)
	$(WGET) $(CMAKE_CUDA_TARBALL_URL)
	@echo '$(CMAKE_CUDA_TARBALL_SHA256)  $(CMAKE_CUDA_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-cmake: .fetch-cmake-$(CMAKE_CUDA_VERSION)

.extract-cmake-$(CMAKE_CUDA_VERSION): .fetch-cmake-$(CMAKE_CUDA_VERSION)
	$(RM) $(CMAKE_SOURCE_DIR)
	$(TAR) -xzf $(CMAKE_CUDA_TARBALL)
	@$(TOUCH) $@

extract-cmake: .extract-cmake-$(CMAKE_CUDA_VERSION)

.configure-cmake-$(CMAKE_CUDA_VERSION): .extract-cmake-$(CMAKE_CUDA_VERSION)
	$(MKDIR) $(CMAKE_BUILD_DIR)
	cd $(CMAKE_BUILD_DIR) && ../configure --prefix=$(CMAKE_INSTALL_DIR)
	@$(TOUCH) $@

configure-cmake: .configure-cmake-$(CMAKE_CUDA_VERSION)

.build-cmake-$(CMAKE_CUDA_VERSION): .configure-cmake-$(CMAKE_CUDA_VERSION)
	cd $(CMAKE_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-cmake: .build-cmake-$(CMAKE_CUDA_VERSION)

install-cmake: .build-cmake-$(CMAKE_CUDA_VERSION)
	cd $(CMAKE_BUILD_DIR) && $(MAKE) install

clean-cmake:
	@$(RM) .build-cmake-$(CMAKE_CUDA_VERSION)
	@$(RM) .configure-cmake-$(CMAKE_CUDA_VERSION)
	$(RM) $(CMAKE_BUILD_DIR)

distclean-cmake: clean-cmake
	@$(RM) .fetch-cmake-$(CMAKE_CUDA_VERSION)
	@$(RM) .extract-cmake-$(CMAKE_CUDA_VERSION)
	$(RM) $(CMAKE_SOURCE_DIR)
	$(RM) $(CMAKE_CUDA_TARBALL)

env-cmake:
	@echo
	@echo '# add CMake $(CMAKE_VERSION) to environment'
	@echo 'export PATH="$(CMAKE_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(CMAKE_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'

##
## Lua
##

ifdef USE_LUA51
LUA_VERSION = 5.1.5
LUA_TARBALL_SHA256 = 2640fc56a795f29d28ef15e13c34a47e223960b0240e8cb0a82d9b0738695333
else
LUA_VERSION = 5.2.2
LUA_TARBALL_SHA256 = 3fd67de3f5ed133bf312906082fa524545c6b9e1b952e8215ffbd27113f49f00
endif
LUA_TARBALL = lua-$(LUA_VERSION).tar.gz
LUA_TARBALL_URL = http://www.lua.org/ftp/$(LUA_TARBALL)
LUA_BUILD_DIR = lua-$(LUA_VERSION)
LUA_INSTALL_DIR = $(PREFIX)/lua-$(LUA_VERSION)
LUA_CFLAGS = -fPIC -O2 -Wall -DLUA_USE_POSIX -DLUA_USE_DLOPEN -DLUA_USE_STRTODHEX -DLUA_USE_AFORMAT -DLUA_USE_LONGLONG
LUA_LIBS = -Wl,-E -ldl -lm

ifdef USE_READLINE
LUA_CFLAGS += -DLUA_USE_READLINE
LUA_LIBS += -lreadline -lncurses
endif

.fetch-lua-$(LUA_VERSION):
	@$(RM) $(LUA_TARBALL)
	$(WGET) $(LUA_TARBALL_URL)
	@echo '$(LUA_TARBALL_SHA256)  $(LUA_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-lua: .fetch-lua-$(LUA_VERSION)

.extract-lua-$(LUA_VERSION): .fetch-lua-$(LUA_VERSION)
	$(RM) $(LUA_BUILD_DIR)
	$(TAR) -xzf $(LUA_TARBALL)
	@$(TOUCH) $@

extract-lua: .extract-lua-$(LUA_VERSION)

.build-lua-$(LUA_VERSION): .extract-lua-$(LUA_VERSION)
	cd $(LUA_BUILD_DIR) && $(MAKE) linux CFLAGS="$(LUA_CFLAGS)" LIBS="$(LUA_LIBS)"
	@$(TOUCH) $@

build-lua: .build-lua-$(LUA_VERSION)

install-lua: .build-lua-$(LUA_VERSION)
	cd $(LUA_BUILD_DIR) && $(MAKE) install INSTALL_TOP=$(LUA_INSTALL_DIR)

clean-lua:
	@$(RM) .build-lua-$(LUA_VERSION)
	@$(RM) .extract-lua-$(LUA_VERSION)
	$(RM) $(LUA_BUILD_DIR)

distclean-lua: clean-lua
	@$(RM) .fetch-lua-$(LUA_VERSION)
	$(RM) $(LUA_TARBALL)

env-lua:
	@echo
	@echo '# add Lua $(LUA_VERSION) to environment'
	@echo 'export PATH="$(LUA_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(LUA_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LUA_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## LuaJIT
##

LUAJIT_VERSION = 2.0.2
LUAJIT_TARBALL = LuaJIT-$(LUAJIT_VERSION).tar.gz
LUAJIT_TARBALL_URL = http://luajit.org/download/$(LUAJIT_TARBALL)
LUAJIT_TARBALL_SHA256 = c05202974a5890e777b181908ac237625b499aece026654d7cc33607e3f46c38
LUAJIT_BUILD_DIR = LuaJIT-$(LUAJIT_VERSION)
LUAJIT_INSTALL_DIR = $(PREFIX)/luajit-$(LUAJIT_VERSION)
LUAJIT_CFLAGS = -fPIC -DLUAJIT_ENABLE_LUA52COMPAT -DLUAJIT_CPU_SSE2

ifdef USE_VALGRIND
LUAJIT_CFLAGS += -DLUAJIT_USE_VALGRIND
endif

.fetch-luajit-$(LUAJIT_VERSION):
	@$(RM) $(LUAJIT_TARBALL)
	$(WGET) $(LUAJIT_TARBALL_URL)
	@echo '$(LUAJIT_TARBALL_SHA256)  $(LUAJIT_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-luajit: .fetch-luajit-$(LUAJIT_VERSION)

.extract-luajit-$(LUAJIT_VERSION): .fetch-luajit-$(LUAJIT_VERSION)
	$(RM) $(LUAJIT_BUILD_DIR)
	$(TAR) -xzf $(LUAJIT_TARBALL)
	@$(TOUCH) $@

extract-luajit: .extract-luajit-$(LUAJIT_VERSION)

.build-luajit-$(LUAJIT_VERSION): .extract-luajit-$(LUAJIT_VERSION)
	cd $(LUAJIT_BUILD_DIR) && $(MAKE) amalg "CFLAGS=$(LUAJIT_CFLAGS)" "PREFIX=$(LUAJIT_INSTALL_DIR)" $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-luajit: .build-luajit-$(LUAJIT_VERSION)

install-luajit: .build-luajit-$(LUAJIT_VERSION)
	cd $(LUAJIT_BUILD_DIR) && $(MAKE) install "PREFIX=$(LUAJIT_INSTALL_DIR)"

clean-luajit:
	@$(RM) .build-luajit-$(LUAJIT_VERSION)
	@$(RM) .extract-luajit-$(LUAJIT_VERSION)
	$(RM) $(LUAJIT_BUILD_DIR)

distclean-luajit: clean-luajit
	@$(RM) .fetch-luajit-$(LUAJIT_VERSION)
	$(RM) $(LUAJIT_TARBALL)

env-luajit:
	@echo
	@echo '# add LuaJIT $(LUAJIT_VERSION) to environment'
	@echo 'export PATH="$(LUAJIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(LUAJIT_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LUAJIT_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## GNU readline wrapper
##

RLWRAP_VERSION = 0.37
RLWRAP_TARBALL = rlwrap-$(RLWRAP_VERSION).tar.gz
RLWRAP_TARBALL_URL = http://utopia.knoware.nl/~hlub/rlwrap/$(RLWRAP_TARBALL)
RLWRAP_TARBALL_SHA256 = 8403a2c184a33ee293a30637afd1362e7dbe0ee642c33b54b2fca68162498bbd
RLWRAP_BUILD_DIR = rlwrap-$(RLWRAP_VERSION)
RLWRAP_INSTALL_DIR = $(PREFIX)/rlwrap-$(RLWRAP_VERSION)

.fetch-rlwrap-$(RLWRAP_VERSION):
	@$(RM) $(RLWRAP_TARBALL)
	$(WGET) $(RLWRAP_TARBALL_URL)
	@echo '$(RLWRAP_TARBALL_SHA256)  $(RLWRAP_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-rlwrap: .fetch-rlwrap-$(RLWRAP_VERSION)

.extract-rlwrap-$(RLWRAP_VERSION): .fetch-rlwrap-$(RLWRAP_VERSION)
	$(RM) $(RLWRAP_BUILD_DIR)
	$(TAR) -xzf $(RLWRAP_TARBALL)
	@$(TOUCH) $@

extract-rlwrap: .extract-rlwrap-$(RLWRAP_VERSION)

.configure-rlwrap-$(RLWRAP_VERSION): .extract-rlwrap-$(RLWRAP_VERSION)
	cd $(RLWRAP_BUILD_DIR) && ./configure -prefix=$(RLWRAP_INSTALL_DIR)
	@$(TOUCH) $@

configure-rlwrap: .configure-rlwrap-$(RLWRAP_VERSION)

.build-rlwrap-$(RLWRAP_VERSION): .configure-rlwrap-$(RLWRAP_VERSION)
	cd $(RLWRAP_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-rlwrap: .build-rlwrap-$(RLWRAP_VERSION)

install-rlwrap: .build-rlwrap-$(RLWRAP_VERSION)
	cd $(RLWRAP_BUILD_DIR) && $(MAKE) install

clean-rlwrap:
	@$(RM) .build-rlwrap-$(RLWRAP_VERSION)
	@$(RM) .configure-rlwrap-$(RLWRAP_VERSION)
	@$(RM) .extract-rlwrap-$(RLWRAP_VERSION)
	$(RM) $(RLWRAP_BUILD_DIR)

distclean-rlwrap: clean-rlwrap
	@$(RM) .fetch-rlwrap-$(RLWRAP_VERSION)
	$(RM) $(RLWRAP_TARBALL)

env-rlwrap:
	@echo
	@echo '# add rlwrap $(RLWRAP_VERSION) to environment'
	@echo 'export PATH="$(RLWRAP_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(RLWRAP_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## luatrace
##

LUATRACE_VERSION = e2a78fa
LUATRACE_GIT_URL = https://github.com/geoffleyland/luatrace.git
LUATRACE_BUILD_DIR = luatrace-$(LUATRACE_VERSION)
LUATRACE_INSTALL_DIR = $(PREFIX)/luatrace-$(LUATRACE_VERSION)

.fetch-luatrace-$(LUATRACE_VERSION):
	@$(RM) $(LUATRACE_BUILD_DIR)
	$(GIT) clone $(LUATRACE_GIT_URL) $(LUATRACE_BUILD_DIR)
	cd $(LUATRACE_BUILD_DIR) && $(GIT) checkout $(LUATRACE_VERSION)
	@$(TOUCH) $@

fetch-luatrace: .fetch-luatrace-$(LUATRACE_VERSION)

install-luatrace: .fetch-luatrace-$(LUATRACE_VERSION)
	cd $(LUATRACE_BUILD_DIR) && $(MAKE) install "LUA=$(LUAJIT_INSTALL_DIR)/bin/luatrace"

clean-luatrace:
	@$(RM) .build-luatrace-$(LUATRACE_VERSION)
	cd $(LUATRACE_BUILD_DIR) && $(MAKE) clean

distclean-luatrace: clean-luatrace
	@$(RM) .fetch-luatrace-$(LUATRACE_VERSION)
	$(RM) $(LUATRACE_BUILD_DIR)

##
## Boost C++ libraries with C++11 ABI
##

BOOST_VERSION = 1.54.0
BOOST_RELEASE = 1_54_0
BOOST_ABI = c++11
BOOST_TARBALL = boost_$(BOOST_RELEASE).tar.bz2
BOOST_TARBALL_URL = http://sourceforge.net/projects/boost/files/boost/$(BOOST_VERSION)/$(BOOST_TARBALL)
BOOST_TARBALL_SHA256 = 047e927de336af106a24bceba30069980c191529fd76b8dff8eb9a328b48ae1d
BOOST_PATCH = boost_$(BOOST_RELEASE).patch
BOOST_PATCH_SHA256 = 28771aa1f84f8bf03360f2ba9189e714025554d4555cb01331bf590f13846db5
BOOST_BUILD_DIR = boost_$(BOOST_RELEASE)
BOOST_INSTALL_DIR = $(PREFIX)/boost_$(BOOST_RELEASE)-$(BOOST_ABI)
BOOST_BUILD_FLAGS = threading=multi variant=release --layout=tagged toolset=gcc cxxflags="-fPIC -std=$(BOOST_ABI)" dll-path=$(BOOST_INSTALL_DIR)/lib

ifndef USE_BZIP2
BOOST_BUILD_FLAGS += -sNO_BZIP2=1
endif
ifndef USE_PYTHON
BOOST_BUILD_FLAGS += --without-python
endif

define BOOST_PATCH_FILE
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
EIH35HBOKor5/gbndk+xdxkAAA==
endef
export BOOST_PATCH_FILE

.fetch-boost-$(BOOST_VERSION):
	@$(RM) $(BOOST_TARBALL)
	@$(RM) $(BOOST_PATCH)
	$(WGET) $(BOOST_TARBALL_URL)
	echo "$$BOOST_PATCH_FILE" | $(BASE64) | $(GUNZIP) > $(BOOST_PATCH)
	@echo '$(BOOST_TARBALL_SHA256)  $(BOOST_TARBALL)' | $(SHA256SUM)
	@echo '$(BOOST_PATCH_SHA256)  $(BOOST_PATCH)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-boost: .fetch-boost-$(BOOST_VERSION)

.extract-boost-$(BOOST_VERSION): .fetch-boost-$(BOOST_VERSION)
	$(TAR) -xjf $(BOOST_TARBALL)
	cd $(BOOST_BUILD_DIR) && $(PATCH) -p0 < ../$(BOOST_PATCH)
	@$(TOUCH) $@

extract-boost: .extract-boost-$(BOOST_VERSION)

.configure-boost-$(BOOST_VERSION): .extract-boost-$(BOOST_VERSION)
	cd $(BOOST_BUILD_DIR) && ./bootstrap.sh
	@$(TOUCH) $@

configure-boost: .configure-boost-$(BOOST_VERSION)

.build-boost-$(BOOST_VERSION): .configure-boost-$(BOOST_VERSION)
	cd $(BOOST_BUILD_DIR) && ./bjam $(BOOST_BUILD_FLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-boost: .build-boost-$(BOOST_VERSION)

install-boost: .build-boost-$(BOOST_VERSION)
	cd $(BOOST_BUILD_DIR) && ./bjam $(BOOST_BUILD_FLAGS) install --prefix=$(BOOST_INSTALL_DIR)

clean-boost:
	$(RM) $(BOOST_BUILD_DIR)
	@$(RM) .build-boost-$(BOOST_VERSION)
	@$(RM) .configure-boost-$(BOOST_VERSION)
	@$(RM) .extract-boost-$(BOOST_VERSION)

distclean-boost: clean-boost
	@$(RM) .fetch-boost-$(BOOST_VERSION)
	$(RM) $(BOOST_TARBALL)
	$(RM) $(BOOST_PATCH)

env-boost:
	@echo
	@echo '# add Boost $(BOOST_VERSION) to environment'
	@echo 'export PATH="$(BOOST_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(BOOST_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export PYTHONPATH="$(BOOST_INSTALL_DIR)/lib$${PYTHONPATH+:$$PYTHONPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(BOOST_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## HDF5 C++ library
##
HDF5_VERSION = 1.8.11
HDF5_TARBALL = hdf5-$(HDF5_VERSION).tar.bz2
HDF5_TARBALL_URL = http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_VERSION)/src/$(HDF5_TARBALL)
HDF5_TARBALL_SHA256 = 5ea1ba82fc77350ee628b795ae6ede05feeaf1c6b437911a9478de456600cafb
HDF5_BUILD_DIR = hdf5-$(HDF5_VERSION)
HDF5_INSTALL_DIR = $(PREFIX)/hdf5-$(HDF5_VERSION)
HDF5_CONFIGURE_FLAGS = --enable-parallel --enable-shared --disable-deprecated-symbols  --enable-cxx
HDF5_CFLAGS = -fPIC
HDF5_CXXFLAGS = -fPIC

.fetch-hdf5-$(HDF5_VERSION):
	@$(RM) $(HDF5_TARBALL)
	$(WGET) $(HDF5_TARBALL_URL)
	@echo '$(HDF5_TARBALL_SHA256)  $(HDF5_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-hdf5: .fetch-hdf5-$(HDF5_VERSION)

.extract-hdf5-$(HDF5_VERSION): .fetch-hdf5-$(HDF5_VERSION)
	$(RM) $(HDF5_BUILD_DIR)
	$(TAR) -xjf $(HDF5_TARBALL)
	@$(TOUCH) $@

extract-hdf5: .extract-hdf5-$(HDF5_VERSION)

.configure-hdf5-$(HDF5_VERSION): .extract-hdf5-$(HDF5_VERSION)
	cd $(HDF5_BUILD_DIR) && CFLAGS="$(HDF5_CFLAGS)" CXXFLAGS="$(HDF5_CXXFLAGS)" ./configure $(HDF5_CONFIGURE_FLAGS) --prefix=$(HDF5_INSTALL_DIR)
	@$(TOUCH) $@

configure-hdf5: .configure-hdf5-$(HDF5_VERSION)

.build-hdf5-$(HDF5_VERSION): .configure-hdf5-$(HDF5_VERSION)
	cd $(HDF5_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-hdf5: .build-hdf5-$(HDF5_VERSION)

install-hdf5: .build-hdf5-$(HDF5_VERSION)
	cd $(HDF5_BUILD_DIR) && $(MAKE) install

clean-hdf5:
	@$(RM) .build-hdf5-$(HDF5_VERSION)
	@$(RM) .configure-hdf5-$(HDF5_VERSION)
	@$(RM) .extract-hdf5-$(HDF5_VERSION)
	$(RM) $(HDF5_BUILD_DIR)

distclean-hdf5: clean-hdf5
	@$(RM) .fetch-hdf5-$(HDF5_VERSION)
	$(RM) $(HDF5_TARBALL)

env-hdf5:
	@echo
	@echo '# add HDF5 $(HDF5_VERSION) to environment'
	@echo 'export PATH="$(HDF5_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(HDF5_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(HDF5_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Git version control
##
GIT_VERSION = 1.8.4
GIT_TARBALL = git-$(GIT_VERSION).tar.gz
GIT_TARBALL_URL = http://git-core.googlecode.com/files/$(GIT_TARBALL)
GIT_TARBALL_SHA256 = 51e8522c256ef8091c6fd5846b9cef8ba4e764819193b5b6cec570f530e6cb94
GIT_MANPAGES_TARBALL = git-manpages-$(GIT_VERSION).tar.gz
GIT_MANPAGES_TARBALL_URL = http://git-core.googlecode.com/files/$(GIT_MANPAGES_TARBALL)
GIT_MANPAGES_TARBALL_SHA256 = 854552e693b22a7ff1c25d2d710d87c8f7dad14888c946c10be60893fdee091f
GIT_BUILD_DIR = git-$(GIT_VERSION)
GIT_CONFIGURE_FLAGS = --without-python
GIT_INSTALL_DIR = $(PREFIX)/git-$(GIT_VERSION)

.fetch-git-$(GIT_VERSION):
	@$(RM) $(GIT_TARBALL)
	@$(RM) $(GIT_MANPAGES_TARBALL)
	$(WGET) $(GIT_TARBALL_URL)
	$(WGET) $(GIT_MANPAGES_TARBALL_URL)
	@echo '$(GIT_TARBALL_SHA256)  $(GIT_TARBALL)' | $(SHA256SUM)
	@echo '$(GIT_MANPAGES_TARBALL_SHA256)  $(GIT_MANPAGES_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-git: .fetch-git-$(GIT_VERSION)

.extract-git-$(GIT_VERSION): .fetch-git-$(GIT_VERSION)
	$(RM) $(GIT_BUILD_DIR)
	$(TAR) -xzf $(GIT_TARBALL)
	@$(TOUCH) $@

extract-git: .extract-git-$(GIT_VERSION)

.configure-git-$(GIT_VERSION): .extract-git-$(GIT_VERSION)
	cd $(GIT_BUILD_DIR) && ./configure $(GIT_CONFIGURE_FLAGS) --prefix=$(GIT_INSTALL_DIR)
	@$(TOUCH) $@

configure-git: .configure-git-$(GIT_VERSION)

.build-git-$(GIT_VERSION): .configure-git-$(GIT_VERSION)
	cd $(GIT_BUILD_DIR) && $(MAKE) PREFIX=$(GIT_INSTALL_DIR) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-git: .build-git-$(GIT_VERSION)

install-git: .build-git-$(GIT_VERSION)
	cd $(GIT_BUILD_DIR) && $(MAKE) PREFIX=$(GIT_INSTALL_DIR) install
	install -d $(GIT_INSTALL_DIR)/share/man
	cd $(GIT_INSTALL_DIR)/share/man && $(TAR) -xzf $(CURDIR)/$(GIT_MANPAGES_TARBALL)

clean-git:
	@$(RM) .build-git-$(GIT_VERSION)
	@$(RM) .configure-git-$(GIT_VERSION)
	@$(RM) .extract-git-$(GIT_VERSION)
	$(RM) $(GIT_BUILD_DIR)

distclean-git: clean-git
	@$(RM) .fetch-git-$(GIT_VERSION)
	$(RM) $(GIT_TARBALL)
	$(RM) $(GIT_MANPAGES_TARBALL)

env-git:
	@echo
	@echo '# add Git $(GIT_VERSION) to environment'
	@echo 'export PATH="$(GIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GIT_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## python-sphinx
##

PYTHON_SPHINX_VERSION = 1.1.3
PYTHON_SPHINX_TARBALL = Sphinx-$(PYTHON_SPHINX_VERSION).tar.gz
PYTHON_SPHINX_TARBALL_URL = http://pypi.python.org/packages/source/S/Sphinx/$(PYTHON_SPHINX_TARBALL)
PYTHON_SPHINX_TARBALL_SHA256 = 34dc95b70a2b07a61b5d61034c34b05f82514aab54ad27adedb49cee911bb8e9
PYTHON_SPHINX_BUILD_DIR = Sphinx-$(PYTHON_SPHINX_VERSION)
PYTHON_SPHINX_INSTALL_DIR = $(PREFIX)/python-sphinx-$(PYTHON_SPHINX_VERSION)
PYTHON_SPHINX_PYTHONPATH = $(PYTHON_SPHINX_INSTALL_DIR)/lib/python

.fetch-python-sphinx-$(PYTHON_SPHINX_VERSION):
	@$(RM) $(PYTHON_SPHINX_TARBALL)
	$(WGET) $(PYTHON_SPHINX_TARBALL_URL)
	@echo '$(PYTHON_SPHINX_TARBALL_SHA256)  $(PYTHON_SPHINX_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-python-sphinx: .fetch-python-sphinx-$(PYTHON_SPHINX_VERSION)

.extract-python-sphinx-$(PYTHON_SPHINX_VERSION): .fetch-python-sphinx-$(PYTHON_SPHINX_VERSION)
	$(RM) $(PYTHON_SPHINX_BUILD_DIR)
	$(TAR) -xzf $(PYTHON_SPHINX_TARBALL)
	@$(TOUCH) $@

extract-python-sphinx: .extract-python-sphinx-$(PYTHON_SPHINX_VERSION)

.build-python-sphinx-$(PYTHON_SPHINX_VERSION): .extract-python-sphinx-$(PYTHON_SPHINX_VERSION)
	cd $(PYTHON_SPHINX_BUILD_DIR) && python setup.py build
	@$(TOUCH) $@

build-python-sphinx: .build-python-sphinx-$(PYTHON_SPHINX_VERSION)

install-python-sphinx: .build-python-sphinx-$(PYTHON_SPHINX_VERSION)
	install -d $(PYTHON_SPHINX_PYTHONPATH)
	cd $(PYTHON_SPHINX_BUILD_DIR) && PYTHONPATH="$(PYTHON_SPHINX_PYTHONPATH)$${PYTHONPATH+:$$PYTHONPATH}" python setup.py install --home=$(PYTHON_SPHINX_INSTALL_DIR)

clean-python-sphinx:
	@$(RM) .build-python-sphinx-$(PYTHON_SPHINX_VERSION)
	@$(RM) .extract-python-sphinx-$(PYTHON_SPHINX_VERSION)
	$(RM) $(PYTHON_SPHINX_BUILD_DIR)

distclean-python-sphinx: clean-python-sphinx
	@$(RM) .fetch-python-sphinx-$(PYTHON_SPHINX_VERSION)
	$(RM) $(PYTHON_SPHINX_TARBALL)

env-python-sphinx:
	@echo
	@echo '# add python-sphinx $(PYTHON_SPHINX_VERSION) to environment'
	@echo 'export PATH="$(PYTHON_SPHINX_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export PYTHONPATH="$(PYTHON_SPHINX_PYTHONPATH)$${PYTHONPATH+:$$PYTHONPATH}"'

##
## Graphviz
##

GRAPHVIZ_VERSION = 2.28.0
GRAPHVIZ_TARBALL = graphviz-$(GRAPHVIZ_VERSION).tar.gz
GRAPHVIZ_TARBALL_URL = http://www.graphviz.org/pub/graphviz/stable/SOURCES/$(GRAPHVIZ_TARBALL)
GRAPHVIZ_TARBALL_SHA256 = d3aa7973c578cae4cc26d9d6498c57ed06680cab9a4e940d0357a3c6527afc76
GRAPHVIZ_BUILD_DIR = graphviz-$(GRAPHVIZ_VERSION)
GRAPHVIZ_CONFIGURE_FLAGS = --with-qt=no --enable-swig=no --enable-python=no
GRAPHVIZ_INSTALL_DIR = $(PREFIX)/graphviz-$(GRAPHVIZ_VERSION)

.fetch-graphviz-$(GRAPHVIZ_VERSION):
	@$(RM) $(GRAPHVIZ_TARBALL)
	$(WGET) $(GRAPHVIZ_TARBALL_URL)
	@echo '$(GRAPHVIZ_TARBALL_SHA256)  $(GRAPHVIZ_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-graphviz: .fetch-graphviz-$(GRAPHVIZ_VERSION)

.extract-graphviz-$(GRAPHVIZ_VERSION): .fetch-graphviz-$(GRAPHVIZ_VERSION)
	$(RM) $(GRAPHVIZ_BUILD_DIR)
	$(TAR) -xzf $(GRAPHVIZ_TARBALL)
	@$(TOUCH) $@

extract-graphviz: .extract-graphviz-$(GRAPHVIZ_VERSION)

.configure-graphviz-$(GRAPHVIZ_VERSION): .extract-graphviz-$(GRAPHVIZ_VERSION)
	cd $(GRAPHVIZ_BUILD_DIR) && ./configure $(GRAPHVIZ_CONFIGURE_FLAGS) --prefix=$(GRAPHVIZ_INSTALL_DIR)
	@$(TOUCH) $@

configure-graphviz: .configure-graphviz-$(GRAPHVIZ_VERSION)

.build-graphviz-$(GRAPHVIZ_VERSION): .configure-graphviz-$(GRAPHVIZ_VERSION)
	cd $(GRAPHVIZ_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-graphviz: .build-graphviz-$(GRAPHVIZ_VERSION)

install-graphviz: .build-graphviz-$(GRAPHVIZ_VERSION)
	cd $(GRAPHVIZ_BUILD_DIR) && $(MAKE) install

clean-graphviz:
	@$(RM) .build-graphviz-$(GRAPHVIZ_VERSION)
	@$(RM) .configure-graphviz-$(GRAPHVIZ_VERSION)
	@$(RM) .extract-graphviz-$(GRAPHVIZ_VERSION)
	$(RM) $(GRAPHVIZ_BUILD_DIR)

distclean-graphviz: clean-graphviz
	@$(RM) .fetch-graphviz-$(GRAPHVIZ_VERSION)
	$(RM) $(GRAPHVIZ_TARBALL)

env-graphviz:
	@echo
	@echo '# add Graphviz $(GRAPHVIZ_VERSION) to environment'
	@echo 'export PATH="$(GRAPHVIZ_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GRAPHVIZ_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## Clang C++ compiler
##
CLANG_VERSION = 3.3
LLVM_TARBALL = llvm-$(CLANG_VERSION).src.tar.gz
LLVM_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(LLVM_TARBALL)
LLVM_TARBALL_SHA256 = 68766b1e70d05a25e2f502e997a3cb3937187a3296595cf6e0977d5cd6727578
CLANG_TARBALL = cfe-$(CLANG_VERSION).src.tar.gz
CLANG_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(CLANG_TARBALL)
CLANG_TARBALL_SHA256 = b1b55de4ab3a57d3e0331a83e0284610191c77d924e3446498d9113d08dfb996
CLANG_BUILD_DIR = llvm-$(CLANG_VERSION).src
CLANG_CONFIGURE_FLAGS = --enable-optimized --enable-bindings=none --enable-shared --with-gcc-toolchain=$(GCC_INSTALL_DIR)
CLANG_BUILD_FLAGS = REQUIRES_RTTI=1
CLANG_INSTALL_DIR = $(PREFIX)/clang-$(CLANG_VERSION)

.fetch-clang-$(CLANG_VERSION):
	@$(RM) $(LLVM_TARBALL)
	@$(RM) $(CLANG_TARBALL)
	$(WGET) $(LLVM_TARBALL_URL)
	$(WGET) $(CLANG_TARBALL_URL)
	@echo '$(LLVM_TARBALL_SHA256)  $(LLVM_TARBALL)' | $(SHA256SUM)
	@echo '$(CLANG_TARBALL_SHA256)  $(CLANG_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-clang: .fetch-clang-$(CLANG_VERSION)

.extract-clang-$(CLANG_VERSION): .fetch-clang-$(CLANG_VERSION)
	$(RM) $(CLANG_BUILD_DIR)
	$(TAR) -xzf $(LLVM_TARBALL)
	cd $(CLANG_BUILD_DIR)/tools && $(TAR) -xzf $(CURDIR)/$(CLANG_TARBALL) && mv cfe-$(CLANG_VERSION).src clang
	@$(TOUCH) $@

extract-clang: .extract-clang-$(CLANG_VERSION)

.configure-clang-$(CLANG_VERSION): .extract-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && ./configure $(CLANG_CONFIGURE_FLAGS) --prefix=$(CLANG_INSTALL_DIR)
	@$(TOUCH) $@

configure-clang: .configure-clang-$(CLANG_VERSION)

.build-clang-$(CLANG_VERSION): .configure-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && $(MAKE) $(CLANG_BUILD_FLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-clang: .build-clang-$(CLANG_VERSION)

install-clang: .build-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && $(MAKE) install

clean-clang:
	@$(RM) .build-clang-$(CLANG_VERSION)
	@$(RM) .configure-clang-$(CLANG_VERSION)
	@$(RM) .extract-clang-$(CLANG_VERSION)
	$(RM) $(CLANG_BUILD_DIR)

distclean-clang: clean-clang
	@$(RM) .fetch-clang-$(CLANG_VERSION)
	$(RM) $(LLVM_TARBALL)
	$(RM) $(CLANG_TARBALL)

env-clang:
	@echo
	@echo '# add Clang $(CLANG_VERSION) to environment'
	@echo 'export PATH="$(CLANG_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(CLANG_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'
	@echo 'export LD_LIBRARY_PATH="$(CLANG_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'

##
## GNU Parallel
##

GNU_PARALLEL_VERSION = 20120622
GNU_PARALLEL_TARBALL = parallel-$(GNU_PARALLEL_VERSION).tar.bz2
GNU_PARALLEL_TARBALL_URL = http://ftp.gnu.org/gnu/parallel/$(GNU_PARALLEL_TARBALL)
GNU_PARALLEL_TARBALL_SHA256 = bfdbb4a66835eacce2bffde23f5e0f12b4e99ced5e348b5f69b7aa97d0123869
GNU_PARALLEL_BUILD_DIR = parallel-$(GNU_PARALLEL_VERSION)
GNU_PARALLEL_INSTALL_DIR = $(PREFIX)/parallel-$(GNU_PARALLEL_VERSION)

.fetch-gnu-parallel-$(GNU_PARALLEL_VERSION):
	@$(RM) $(GNU_PARALLEL_TARBALL)
	$(WGET) $(GNU_PARALLEL_TARBALL_URL)
	@echo '$(GNU_PARALLEL_TARBALL_SHA256)  $(GNU_PARALLEL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gnu-parallel: .fetch-gnu-parallel-$(GNU_PARALLEL_VERSION)

.extract-gnu-parallel-$(GNU_PARALLEL_VERSION): .fetch-gnu-parallel-$(GNU_PARALLEL_VERSION)
	$(RM) $(GNU_PARALLEL_BUILD_DIR)
	$(TAR) -xjf $(GNU_PARALLEL_TARBALL)
	@$(TOUCH) $@

extract-gnu-parallel: .extract-gnu-parallel-$(GNU_PARALLEL_VERSION)

.configure-gnu-parallel-$(GNU_PARALLEL_VERSION): .extract-gnu-parallel-$(GNU_PARALLEL_VERSION)
	cd $(GNU_PARALLEL_BUILD_DIR) && ./configure --prefix=$(GNU_PARALLEL_INSTALL_DIR)
	@$(TOUCH) $@

configure-gnu-parallel: .configure-gnu-parallel-$(GNU_PARALLEL_VERSION)

.build-gnu-parallel-$(GNU_PARALLEL_VERSION): .configure-gnu-parallel-$(GNU_PARALLEL_VERSION)
	cd $(GNU_PARALLEL_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gnu-parallel: .build-gnu-parallel-$(GNU_PARALLEL_VERSION)

install-gnu-parallel: .build-gnu-parallel-$(GNU_PARALLEL_VERSION)
	cd $(GNU_PARALLEL_BUILD_DIR) && $(MAKE) install

clean-gnu-parallel:
	@$(RM) .build-gnu-parallel-$(GNU_PARALLEL_VERSION)
	@$(RM) .configure-gnu-parallel-$(GNU_PARALLEL_VERSION)
	@$(RM) .extract-gnu-parallel-$(GNU_PARALLEL_VERSION)
	$(RM) $(GNU_PARALLEL_BUILD_DIR)

distclean-gnu-parallel: clean-gnu-parallel
	@$(RM) .fetch-gnu-parallel-$(GNU_PARALLEL_VERSION)
	$(RM) $(GNU_PARALLEL_TARBALL)

env-gnu-parallel:
	@echo
	@echo '# add GNU Parallel $(GNU_PARALLEL_VERSION) to environment'
	@echo 'export PATH="$(GNU_PARALLEL_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GNU_PARALLEL_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## GMP (GNU Multiple Precision Arithmetic Library)
##

GMP_VERSION = 5.0.4
GMP_TARBALL = gmp-$(GMP_VERSION).tar.bz2
GMP_TARBALL_URL = http://ftp.gnu.org/gnu/gmp/$(GMP_TARBALL)
GMP_TARBALL_SHA256 = 35d4aade3e4bdf0915c944599b10d23f108ffedf6c3188aeec52221c5cf9a06f
GMP_BUILD_DIR = gmp-$(GMP_VERSION)
GMP_INSTALL_DIR = $(CURDIR)/.gmp-$(GMP_VERSION)

.fetch-gmp-$(GMP_VERSION):
	@$(RM) $(GMP_TARBALL)
	$(WGET) $(GMP_TARBALL_URL)
	@echo '$(GMP_TARBALL_SHA256)  $(GMP_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gmp: .fetch-gmp-$(GMP_VERSION)

.extract-gmp-$(GMP_VERSION): .fetch-gmp-$(GMP_VERSION)
	$(RM) $(GMP_BUILD_DIR)
	$(TAR) -xjf $(GMP_TARBALL)
	@$(TOUCH) $@

extract-gmp: .extract-gmp-$(GMP_VERSION)

.configure-gmp-$(GMP_VERSION): .extract-gmp-$(GMP_VERSION)
	cd $(GMP_BUILD_DIR) && CFLAGS=-fPIC ./configure --prefix=$(GMP_INSTALL_DIR) --libdir=$(GMP_INSTALL_DIR)/lib --enable-cxx --disable-shared
	@$(TOUCH) $@

configure-gmp: .configure-gmp-$(GMP_VERSION)

.build-gmp-$(GMP_VERSION): .configure-gmp-$(GMP_VERSION)
	cd $(GMP_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gmp: .build-gmp-$(GMP_VERSION)

.install-gmp-$(GMP_VERSION): .build-gmp-$(GMP_VERSION)
	cd $(GMP_BUILD_DIR) && $(MAKE) install
	@$(TOUCH) $@

clean-gmp:
	@$(RM) .build-gmp-$(GMP_VERSION)
	@$(RM) .configure-gmp-$(GMP_VERSION)
	@$(RM) .extract-gmp-$(GMP_VERSION)
	$(RM) $(GMP_BUILD_DIR)

distclean-gmp: clean-gmp
	@$(RM) .fetch-gmp-$(GMP_VERSION)
	$(RM) $(GMP_TARBALL)

##
## MPFR (Multiple-precision floating-point computations with correct rounding)
##

MPFR_VERSION = 3.1.0
MPFR_TARBALL = mpfr-$(MPFR_VERSION).tar.bz2
MPFR_TARBALL_URL = http://www.mpfr.org/mpfr-$(MPFR_VERSION)/$(MPFR_TARBALL)
MPFR_TARBALL_SHA256 = 74a7bbbad168dd1cc414f1c9210b8fc16ccfc8e422d34b3371a8978e31eab680
MPFR_BUILD_DIR = mpfr-$(MPFR_VERSION)
MPFR_INSTALL_DIR = $(CURDIR)/.mpfr-$(MPFR_VERSION)

.fetch-mpfr-$(MPFR_VERSION):
	@$(RM) $(MPFR_TARBALL)
	$(WGET) $(MPFR_TARBALL_URL)
	@echo '$(MPFR_TARBALL_SHA256)  $(MPFR_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-mpfr: .fetch-mpfr-$(MPFR_VERSION)

.extract-mpfr-$(MPFR_VERSION): .fetch-mpfr-$(MPFR_VERSION)
	$(RM) $(MPFR_BUILD_DIR)
	$(TAR) -xjf $(MPFR_TARBALL)
	@$(TOUCH) $@

extract-mpfr: .extract-mpfr-$(MPFR_VERSION)

.configure-mpfr-$(MPFR_VERSION): .extract-mpfr-$(MPFR_VERSION) .install-gmp-$(GMP_VERSION)
	cd $(MPFR_BUILD_DIR) && ./configure --prefix=$(MPFR_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --disable-shared
	@$(TOUCH) $@

configure-mpfr: .configure-mpfr-$(MPFR_VERSION)

.build-mpfr-$(MPFR_VERSION): .configure-mpfr-$(MPFR_VERSION)
	cd $(MPFR_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-mpfr: .build-mpfr-$(MPFR_VERSION)

.install-mpfr-$(MPFR_VERSION): .build-mpfr-$(MPFR_VERSION)
	cd $(MPFR_BUILD_DIR) && $(MAKE) install
	@$(TOUCH) $@

clean-mpfr:
	@$(RM) .build-mpfr-$(MPFR_VERSION)
	@$(RM) .configure-mpfr-$(MPFR_VERSION)
	@$(RM) .extract-mpfr-$(MPFR_VERSION)
	$(RM) $(MPFR_BUILD_DIR)

distclean-mpfr: clean-mpfr
	@$(RM) .fetch-mpfr-$(MPFR_VERSION)
	$(RM) $(MPFR_TARBALL)

##
## MPC (arithmetic of complex numbers with arbitrarily high precision and correct rounding)
##

MPC_VERSION = 0.9
MPC_TARBALL = mpc-$(MPC_VERSION).tar.gz
MPC_TARBALL_URL = http://www.multiprecision.org/mpc/download/$(MPC_TARBALL)
MPC_TARBALL_SHA256 = fd3efe422f0d454592059e80f2c00d1a2e381bf2beda424c5094abd4deb049ac
MPC_BUILD_DIR = mpc-$(MPC_VERSION)
MPC_INSTALL_DIR = $(CURDIR)/.mpc-$(MPC_VERSION)

.fetch-mpc-$(MPC_VERSION):
	@$(RM) $(MPC_TARBALL)
	$(WGET) $(MPC_TARBALL_URL)
	@echo '$(MPC_TARBALL_SHA256)  $(MPC_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-mpc: .fetch-mpc-$(MPC_VERSION)

.extract-mpc-$(MPC_VERSION): .fetch-mpc-$(MPC_VERSION)
	$(RM) $(MPC_BUILD_DIR)
	$(TAR) -xzf $(MPC_TARBALL)
	@$(TOUCH) $@

extract-mpc: .extract-mpc-$(MPC_VERSION)

.configure-mpc-$(MPC_VERSION): .extract-mpc-$(MPC_VERSION) .install-gmp-$(GMP_VERSION) .install-mpfr-$(MPFR_VERSION)
	cd $(MPC_BUILD_DIR) && ./configure --prefix=$(MPC_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --with-mpfr=$(MPFR_INSTALL_DIR) --disable-shared
	@$(TOUCH) $@

configure-mpc: .configure-mpc-$(MPC_VERSION)

.build-mpc-$(MPC_VERSION): .configure-mpc-$(MPC_VERSION)
	cd $(MPC_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-mpc: .build-mpc-$(MPC_VERSION)

.install-mpc-$(MPC_VERSION): .build-mpc-$(MPC_VERSION)
	cd $(MPC_BUILD_DIR) && $(MAKE) install
	@$(TOUCH) $@

clean-mpc:
	@$(RM) .build-mpc-$(MPC_VERSION)
	@$(RM) .configure-mpc-$(MPC_VERSION)
	@$(RM) .extract-mpc-$(MPC_VERSION)
	$(RM) $(MPC_BUILD_DIR)

distclean-mpc: clean-mpc
	@$(RM) .fetch-mpc-$(MPC_VERSION)
	$(RM) $(MPC_TARBALL)

##
## PPL (Parma Polyhedra Library)
##

PPL_VERSION = 0.11.2
PPL_TARBALL = ppl-$(PPL_VERSION).tar.bz2
PPL_TARBALL_URL = ftp://ftp.cs.unipr.it/pub/ppl/releases/$(PPL_VERSION)/$(PPL_TARBALL)
PPL_TARBALL_SHA256 = e3fbd1c19ef44c6f020951807cdb6fc6a8153cd3a5c53b0ab9cf4c4f6e8cbbeb
PPL_BUILD_DIR = ppl-$(PPL_VERSION)
PPL_INSTALL_DIR = $(CURDIR)/.ppl-$(PPL_VERSION)

.fetch-ppl-$(PPL_VERSION):
	@$(RM) $(PPL_TARBALL)
	$(WGET) $(PPL_TARBALL_URL)
	@echo '$(PPL_TARBALL_SHA256)  $(PPL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-ppl: .fetch-ppl-$(PPL_VERSION)

.extract-ppl-$(PPL_VERSION): .fetch-ppl-$(PPL_VERSION)
	$(RM) $(PPL_BUILD_DIR)
	$(TAR) -xjf $(PPL_TARBALL)
	@$(TOUCH) $@

extract-ppl: .extract-ppl-$(PPL_VERSION)

.configure-ppl-$(PPL_VERSION): .extract-ppl-$(PPL_VERSION) .install-gmp-$(GMP_VERSION)
	cd $(PPL_BUILD_DIR) && ./configure --prefix=$(PPL_INSTALL_DIR) --libdir=$(PPL_INSTALL_DIR)/lib --with-gmp=$(GMP_INSTALL_DIR) --with-gmp-library=$(GMP_INSTALL_DIR)/lib --enable-interfaces=c,cxx --disable-shared --disable-watchdog
	@$(TOUCH) $@

configure-ppl: .configure-ppl-$(PPL_VERSION)

.build-ppl-$(PPL_VERSION): .configure-ppl-$(PPL_VERSION)
	cd $(PPL_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-ppl: .build-ppl-$(PPL_VERSION)

.install-ppl-$(PPL_VERSION): .build-ppl-$(PPL_VERSION)
	cd $(PPL_BUILD_DIR) && $(MAKE) install
	@$(TOUCH) $@

clean-ppl:
	@$(RM) .build-ppl-$(PPL_VERSION)
	@$(RM) .configure-ppl-$(PPL_VERSION)
	@$(RM) .extract-ppl-$(PPL_VERSION)
	$(RM) $(PPL_BUILD_DIR)

distclean-ppl: clean-ppl
	@$(RM) .fetch-ppl-$(PPL_VERSION)
	$(RM) $(PPL_TARBALL)

##
## CLooG-PPL
##

CLOOG_PPL_VERSION = 0.15.11
CLOOG_PPL_TARBALL = cloog-ppl-$(CLOOG_PPL_VERSION).tar.gz
CLOOG_PPL_TARBALL_URL = ftp://ftp.fu-berlin.de/unix/languages/gcc/infrastructure/$(CLOOG_PPL_TARBALL)
CLOOG_PPL_TARBALL_SHA256 = 7cd634d0b2b401b04096b545915ac67f883556e9a524e8e803a6bf6217a84d5f
CLOOG_PPL_BUILD_DIR = cloog-ppl-$(CLOOG_PPL_VERSION)
CLOOG_PPL_INSTALL_DIR = $(CURDIR)/.cloog-ppl-$(CLOOG_PPL_VERSION)

.fetch-cloog-ppl-$(CLOOG_PPL_VERSION):
	@$(RM) $(CLOOG_PPL_TARBALL)
	$(WGET) $(CLOOG_PPL_TARBALL_URL)
	@echo '$(CLOOG_PPL_TARBALL_SHA256)  $(CLOOG_PPL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-cloog-ppl: .fetch-cloog-ppl-$(CLOOG_PPL_VERSION)

.extract-cloog-ppl-$(CLOOG_PPL_VERSION): .fetch-cloog-ppl-$(CLOOG_PPL_VERSION)
	$(RM) $(CLOOG_PPL_BUILD_DIR)
	$(TAR) -xzf $(CLOOG_PPL_TARBALL)
	@$(TOUCH) $@

extract-cloog-ppl: .extract-cloog-ppl-$(CLOOG_PPL_VERSION)

.configure-cloog-ppl-$(CLOOG_PPL_VERSION): .extract-cloog-ppl-$(CLOOG_PPL_VERSION) .install-ppl-$(PPL_VERSION) .install-gmp-$(GMP_VERSION)
	cd $(CLOOG_PPL_BUILD_DIR) && LDFLAGS="-L$(PPL_INSTALL_DIR)/lib -lm" ./configure --prefix=$(CLOOG_PPL_INSTALL_DIR) --libdir=$(CLOOG_PPL_INSTALL_DIR)/lib --with-ppl=$(PPL_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --with-gmp-library=$(GMP_INSTALL_DIR)/lib --disable-shared --with-host-libstdcxx=-lstdc++
	@$(TOUCH) $@

configure-cloog-ppl: .configure-cloog-ppl-$(CLOOG_PPL_VERSION)

.build-cloog-ppl-$(CLOOG_PPL_VERSION): .configure-cloog-ppl-$(CLOOG_PPL_VERSION)
	cd $(CLOOG_PPL_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-cloog-ppl: .build-cloog-ppl-$(CLOOG_PPL_VERSION)

.install-cloog-ppl-$(CLOOG_PPL_VERSION): .build-cloog-ppl-$(CLOOG_PPL_VERSION)
	cd $(CLOOG_PPL_BUILD_DIR) && $(MAKE) install
	@$(TOUCH) $@

clean-cloog-ppl:
	@$(RM) .build-cloog-ppl-$(CLOOG_PPL_VERSION)
	@$(RM) .configure-cloog-ppl-$(CLOOG_PPL_VERSION)
	@$(RM) .extract-cloog-ppl-$(CLOOG_PPL_VERSION)
	$(RM) $(CLOOG_PPL_BUILD_DIR)

distclean-cloog-ppl: clean-cloog-ppl
	@$(RM) .fetch-cloog-ppl-$(CLOOG_PPL_VERSION)
	$(RM) $(CLOOG_PPL_TARBALL)

##
## GCC (GNU Compiler Collection)
##

GCC_VERSION = 4.8.2
GCC_TARBALL = gcc-$(GCC_VERSION).tar.bz2
GCC_TARBALL_URL = http://ftp.gwdg.de/pub/misc/gcc/releases/gcc-$(GCC_VERSION)/$(GCC_TARBALL)
GCC_TARBALL_SHA256 = 09dc2276c73424bbbfda1dbddc62bbbf900c9f185acf7f3e1d773ce2d7e3cdc8
GCC_BUILD_DIR = gcc-$(GCC_VERSION)
GCC_BUILD_FLAGS = --enable-cxx-flags=-fPIC --enable-languages=c,c++,fortran,lto --disable-multilib
GCC_INSTALL_DIR = $(PREFIX)/gcc-$(GCC_VERSION)

.fetch-gcc-$(GCC_VERSION):
	@$(RM) $(GCC_TARBALL)
	$(WGET) $(GCC_TARBALL_URL)
	@echo '$(GCC_TARBALL_SHA256)  $(GCC_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gcc: .fetch-gcc-$(GCC_VERSION) .fetch-gmp-$(GMP_VERSION) .fetch-mpfr-$(MPFR_VERSION) .fetch-mpc-$(MPC_VERSION) .fetch-ppl-$(PPL_VERSION) .fetch-cloog-ppl-$(CLOOG_PPL_VERSION)

.extract-gcc-$(GCC_VERSION): .fetch-gcc-$(GCC_VERSION)
	$(RM) $(GCC_BUILD_DIR)
	$(TAR) -xjf $(GCC_TARBALL)
	@$(TOUCH) $@

extract-gcc: .extract-gcc-$(GCC_VERSION)

.configure-gcc-$(GCC_VERSION): .extract-gcc-$(GCC_VERSION) .install-gmp-$(GMP_VERSION) .install-mpfr-$(MPFR_VERSION) .install-mpc-$(MPC_VERSION) .install-ppl-$(PPL_VERSION) .install-cloog-ppl-$(CLOOG_PPL_VERSION)
	cd $(GCC_BUILD_DIR) && ac_cv_lib_pwl_PWL_handle_timeout=no ./configure --prefix=$(GCC_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --with-gmp-library=$(GMP_INSTALL_DIR)/lib --with-mpfr=$(MPFR_INSTALL_DIR) --with-mpc=$(MPC_INSTALL_DIR) --with-ppl=$(PPL_INSTALL_DIR) --with-ppl-library=$(PPL_INSTALL_DIR)/lib --with-cloog=$(CLOOG_PPL_INSTALL_DIR) --enable-build-with-cxx $(GCC_BUILD_FLAGS)
	@$(TOUCH) $@

configure-gcc: .configure-gcc-$(GCC_VERSION)

.build-gcc-$(GCC_VERSION): .configure-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gcc: .build-gcc-$(GCC_VERSION)

install-gcc: .build-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && $(MAKE) install
	cd $(GMP_INSTALL_DIR) && $(CP) include `$(GCC_INSTALL_DIR)/bin/gcc -print-file-name=plugin`

clean-gcc:
	@$(RM) .build-gcc-$(GCC_VERSION)
	@$(RM) .configure-gcc-$(GCC_VERSION)
	@$(RM) .extract-gcc-$(GCC_VERSION)
	$(RM) $(GCC_BUILD_DIR)

distclean-gcc: clean-gcc
	@$(RM) .fetch-gcc-$(GCC_VERSION)
	$(RM) $(GCC_TARBALL)

env-gcc:
	@echo
	@echo '# add gcc $(GCC_VERSION) to environment'
	@echo 'export PATH="$(GCC_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GCC_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'
	@echo 'export LD_LIBRARY_PATH="$(GCC_INSTALL_DIR)/lib64$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'

##
## HALMD Highly Accerelated Large-scale Molecular Dynamics
##

HALMD_VERSION = 0.2.1
HALMD_GIT_URL = http://git.halmd.org/halmd.git
HALMD_SOURCE_DIR = halmd-$(HALMD_VERSION)
HALMD_BUILD_DIR = $(HALMD_SOURCE_DIR)/build/release
HALMD_INSTALL_DIR = $(PREFIX)/halmd-$(HALMD_VERSION)
HALMD_BUILD_ENV = 

.fetch-halmd-$(HALMD_VERSION):
	$(RM) $(HALMD_SOURCE_DIR)
	$(GIT) clone $(HALMD_GIT_URL) $(HALMD_SOURCE_DIR)
	cd $(HALMD_SOURCE_DIR) && $(GIT) checkout $(HALMD_VERSION)
	cd $(HALMD_SOURCE_DIR) && $(GIT) submodule update --init
	@$(TOUCH) $@

fetch-halmd: .fetch-halmd-$(HALMD_VERSION)

.configure-halmd-$(HALMD_VERSION): .fetch-halmd-$(HALMD_VERSION)
	$(MKDIR) $(HALMD_BUILD_DIR)
	cd $(HALMD_BUILD_DIR) && $(HALMD_BUILD_ENV) cmake -DCMAKE_INSTALL_PREFIX=$(HALMD_INSTALL_DIR) $(CURDIR)/$(HALMD_SOURCE_DIR)
	@$(TOUCH) $@

configure-halmd: .configure-halmd-$(HALMD_VERSION)

.build-halmd-$(HALMD_VERSION): .configure-halmd-$(HALMD_VERSION)
	cd $(HALMD_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-halmd: .build-halmd-$(HALMD_VERSION)

install-halmd: .build-halmd-$(HALMD_VERSION)
	cd $(HALMD_BUILD_DIR) && $(MAKE) install

clean-halmd:
	@$(RM) .build-halmd-$(HALMD_VERSION)
	@$(RM) .configure-halmd-$(HALMD_VERSION)
	$(RM) $(HALMD_BUILD_DIR)

distclean-halmd: clean-halmd
	@$(RM) .fetch-halmd-$(HALMD_VERSION)
	$(RM) $(HALMD_SOURCE_DIR)

env-halmd:
	@echo
	@echo '# add HALMD $(HALMD_VERSION) to environment'
	@echo 'export PATH="$(HALMD_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(HALMD_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## nvCUDA-tools - A collection of tools for NVIDIA CUDA compute devices
##

NVCUDA_TOOLS_VERSION = master
NVCUDA_TOOLS_GIT_URL = http://github.com/fhoefling/nvcuda-tools.git
NVCUDA_TOOLS_SOURCE_DIR = nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
NVCUDA_TOOLS_BUILD_DIR = $(NVCUDA_TOOLS_SOURCE_DIR)/build/release
NVCUDA_TOOLS_INSTALL_DIR = $(PREFIX)/nvcuda-tools-$(NVCUDA_TOOLS_VERSION)

.fetch-nvcuda-tools-$(NVCUDA_TOOLS_VERSION):
	$(RM) $(NVCUDA_TOOLS_SOURCE_DIR)
	$(GIT) clone $(NVCUDA_TOOLS_GIT_URL) $(NVCUDA_TOOLS_SOURCE_DIR)
	cd $(NVCUDA_TOOLS_SOURCE_DIR) && $(GIT) checkout $(NVCUDA_TOOLS_VERSION)
	@$(TOUCH) $@

fetch-nvcuda-tools: .fetch-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)

.configure-nvcuda-tools-$(NVCUDA_TOOLS_VERSION): .fetch-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	$(MKDIR) $(NVCUDA_TOOLS_BUILD_DIR)
	cd $(NVCUDA_TOOLS_BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(NVCUDA_TOOLS_INSTALL_DIR) $(CURDIR)/$(NVCUDA_TOOLS_SOURCE_DIR)
	@$(TOUCH) $@

configure-nvcuda-tools: .configure-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)

.build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION): .configure-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	cd $(NVCUDA_TOOLS_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-nvcuda-tools: .build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)

install-nvcuda-tools: .build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	cd $(NVCUDA_TOOLS_BUILD_DIR) && $(MAKE) install

clean-nvcuda-tools:
	@$(RM) .build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	@$(RM) .configure-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	$(RM) $(NVCUDA_TOOLS_BUILD_DIR)

distclean-nvcuda-tools: clean-nvcuda-tools
	@$(RM) .fetch-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	$(RM) $(NVCUDA_TOOLS_SOURCE_DIR)

env-nvcuda-tools:
	@echo
	@echo '# add nvCUDA-tools $(NVCUDA_TOOLS_VERSION) to environment'
	@echo 'export PATH="$(NVCUDA_TOOLS_INSTALL_DIR)/bin$${PATH+:$$PATH}"'

