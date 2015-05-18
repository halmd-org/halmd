#
# Copyright © 2011-2014 Peter Colberg
# Copyright © 2013      Nicolas Höft
# Copyright © 2013-2014 Felix Höfling
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
CTEST     = ctest
GIT       = git
GUNZIP    = gzip -d
MKDIR     = mkdir -p
OPENSSL   = openssl
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

build: build-cmake build-lua build-luajit build-boost build-hdf5 build-git build-python-sphinx build-graphviz build-gdb build-clang build-gcc build-halmd build-nvcuda-tools build-ninja

fetch: fetch-cmake fetch-lua fetch-luajit fetch-boost fetch-hdf5 fetch-git fetch-python-sphinx fetch-graphviz fetch-gdb fetch-clang fetch-gcc fetch-halmd fetch-nvcuda-tools fetch-ninja

install: install-cmake install-lua install-luajit install-boost install-hdf5 install-git install-python-sphinx install-graphviz install-gdb install-clang install-gcc install-halmd install-nvcuda-tools install-ninja

clean: clean-cmake clean-lua clean-luajit clean-boost clean-hdf5 clean-git clean-python-sphinx clean-graphviz clean-gdb clean-clang clean-gcc clean-halmd clean-nvcuda-tools clean-ninja

distclean: distclean-cmake distclean-lua distclean-luajit distclean-boost distclean-hdf5 distclean-git distclean-python-sphinx distclean-graphviz distclean-gdb distclean-clang distclean-gcc distclean-halmd distclean-nvcuda-tools distclean-ninja

env: env-cmake env-lua env-luajit env-boost env-hdf5 env-git env-python-sphinx env-graphviz env-gdb env-clang env-gcc env-halmd env-nvcuda-tools env-ninja

##
## CMake with CMake-CUDA patch
##

CMAKE_CUDA_VERSION = 2.8.12.2-cuda
CMAKE_CUDA_TARBALL = $(CMAKE_CUDA_VERSION).tar.gz
CMAKE_CUDA_TARBALL_SHA256 = e884782cfc9fcb2c3720c51d1924d2a0dacebb10d09a6524e773ada1ba76c093
CMAKE_CUDA_TARBALL_URL = http://git.halmd.org/cmake-cuda/archive/$(CMAKE_CUDA_TARBALL)
CMAKE_SOURCE_DIR = cmake-cuda-$(CMAKE_CUDA_VERSION)
CMAKE_BUILD_DIR = $(CMAKE_SOURCE_DIR)/build
CMAKE_INSTALL_DIR = $(PREFIX)/cmake-$(CMAKE_CUDA_VERSION)

.fetch-cmake-$(CMAKE_CUDA_VERSION):
	@$(RM) $(CMAKE_CUDA_TARBALL)
	$(WGET) $(CMAKE_CUDA_TARBALL_URL) -O $(CMAKE_CUDA_TARBALL)
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
	@echo 'export PATH="$(CMAKE_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(CMAKE_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'

##
## Lua
##

ifdef USE_LUA51
LUA_VERSION = 5.1.5
LUA_TARBALL_SHA256 = 2640fc56a795f29d28ef15e13c34a47e223960b0240e8cb0a82d9b0738695333
else
LUA_VERSION = 5.3.0
LUA_TARBALL_SHA256 = ae4a5eb2d660515eb191bfe3e061f2b8ffe94dce73d32cfd0de090ddcc0ddb01
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
	@echo 'export PATH="$(LUA_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(LUA_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LUA_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## LuaJIT
##

LUAJIT_VERSION = 2.0.4
LUAJIT_TARBALL = LuaJIT-$(LUAJIT_VERSION).tar.gz
LUAJIT_TARBALL_URL = http://luajit.org/download/$(LUAJIT_TARBALL)
LUAJIT_TARBALL_SHA256 = 620fa4eb12375021bef6e4f237cbd2dd5d49e56beb414bee052c746beef1807d
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
	@echo 'export PATH="$(LUAJIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(LUAJIT_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LUAJIT_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## luatrace
##

LUATRACE_VERSION = 6150cfd
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

BOOST_VERSION = 1.58.0
BOOST_RELEASE = 1_58_0
BOOST_ABI = c++11
BOOST_TOOLSET = gcc
BOOST_TARBALL = boost_$(BOOST_RELEASE).tar.bz2
BOOST_TARBALL_URL = http://sourceforge.net/projects/boost/files/boost/$(BOOST_VERSION)/$(BOOST_TARBALL)
BOOST_TARBALL_SHA256 = fdfc204fc33ec79c99b9a74944c3e54bd78be4f7f15e260c0e2700a36dc7d3e5
BOOST_BUILD_DIR = boost_$(BOOST_RELEASE)
BOOST_INSTALL_DIR = $(PREFIX)/boost_$(BOOST_RELEASE)-$(BOOST_ABI)
BOOST_BUILD_FLAGS = threading=multi variant=release --layout=tagged toolset=$(BOOST_TOOLSET) cxxflags="-fPIC -std=$(BOOST_ABI)" dll-path=$(BOOST_INSTALL_DIR)/lib

ifndef USE_BZIP2
BOOST_BUILD_FLAGS += -sNO_BZIP2=1
endif
ifndef USE_PYTHON
BOOST_BUILD_FLAGS += --without-python
endif

.fetch-boost-$(BOOST_VERSION):
	@$(RM) $(BOOST_TARBALL)
	$(WGET) $(BOOST_TARBALL_URL)
	@echo '$(BOOST_TARBALL_SHA256)  $(BOOST_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-boost: .fetch-boost-$(BOOST_VERSION)

.extract-boost-$(BOOST_VERSION): .fetch-boost-$(BOOST_VERSION)
	$(TAR) -xjf $(BOOST_TARBALL)
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

env-boost:
	@echo 'export PATH="$(BOOST_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(BOOST_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export PYTHONPATH="$(BOOST_INSTALL_DIR)/lib$${PYTHONPATH+:$$PYTHONPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(BOOST_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## HDF5 library
##
HDF5_VERSION = 1.8.15
HDF5_TARBALL = hdf5-$(HDF5_VERSION).tar.bz2
HDF5_TARBALL_URL = http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_VERSION)/src/$(HDF5_TARBALL)
HDF5_TARBALL_SHA256 = 0003ac76314b9c3dd819f595eb33890f02f6bf000414013a89f9907ac09e3920
HDF5_BUILD_DIR = hdf5-$(HDF5_VERSION)
HDF5_INSTALL_DIR = $(PREFIX)/hdf5-$(HDF5_VERSION)
HDF5_CONFIGURE_FLAGS = --enable-shared --enable-cxx --enable-fortran
HDF5_CFLAGS = -fPIC
HDF5_CXXFLAGS = -fPIC

ifdef ENABLE_PARALLEL_HDF5
    HDF5_CONFIGURE_FLAGS += --enable-parallel
endif

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
	@echo 'export PATH="$(HDF5_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(HDF5_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(HDF5_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Curl library, needed for Git to use the https protocol
##

CURL_VERSION = 7.38.0
CURL_TARBALL = curl-$(CURL_VERSION).tar.gz
CURL_TARBALL_URL = http://curl.haxx.se/download/$(CURL_TARBALL)
CURL_TARBALL_SHA256 = 5661028aa6532882fa228cd23c99ddbb8b87643dbb1a7ea55c068d34a943dff1
CURL_BUILD_DIR = curl-$(CURL_VERSION)
CURL_INSTALL_DIR = $(CURDIR)/.curl-$(CURL_VERSION)
CURL_CONFIGURE_FLAGS = --with-ssl

.fetch-curl-$(CURL_VERSION):
	@$(RM) $(CURL_TARBALL)
	$(WGET) $(CURL_TARBALL_URL)
	@echo '$(CURL_TARBALL_SHA256)  $(CURL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-curl: .fetch-curl-$(CURL_VERSION)

.extract-curl-$(CURL_VERSION): .fetch-curl-$(CURL_VERSION)
	$(RM) $(CURL_BUILD_DIR)
	$(TAR) -xzf $(CURL_TARBALL)
	@$(TOUCH) $@

extract-curl: .extract-curl-$(CURL_VERSION)

.configure-curl-$(CURL_VERSION): .extract-curl-$(CURL_VERSION)
	cd $(CURL_BUILD_DIR) && ./configure --prefix=$(CURL_INSTALL_DIR) --libdir=$(CURL_INSTALL_DIR)/lib $(CURL_CONFIGURE_FLAGS)
	@$(TOUCH) $@

configure-curl: .configure-curl-$(CURL_VERSION)

.build-curl-$(CURL_VERSION): .configure-curl-$(CURL_VERSION)
	cd $(CURL_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-curl: .build-curl-$(CURL_VERSION)

.install-curl-$(CURL_VERSION): .build-curl-$(CURL_VERSION)
	cd $(CURL_BUILD_DIR) && $(MAKE) install
	@$(TOUCH) $@

clean-curl:
	@$(RM) .build-curl-$(CURL_VERSION)
	@$(RM) .configure-curl-$(CURL_VERSION)
	@$(RM) .extract-curl-$(CURL_VERSION)
	$(RM) $(CURL_BUILD_DIR)

distclean-curl: clean-curl
	@$(RM) .fetch-curl-$(CURL_VERSION)
	@$(RM) .install-curl-$(CURL_VERSION)
	$(RM) $(CURL_TARBALL)
	$(RM) $(CURL_INSTALL_DIR)

##
## Git version control
##
GIT_VERSION = 2.4.0
GIT_TARBALL = git-$(GIT_VERSION).tar.gz
GIT_TARBALL_URL = https://www.kernel.org/pub/software/scm/git/$(GIT_TARBALL)
GIT_TARBALL_SHA256 = d58c766a80d86a5e1846c04c74618e98bfec158734fa1a8904c64740b72b511a
GIT_MANPAGES_TARBALL = git-manpages-$(GIT_VERSION).tar.gz
GIT_MANPAGES_TARBALL_URL = https://www.kernel.org/pub/software/scm/git/$(GIT_MANPAGES_TARBALL)
GIT_MANPAGES_TARBALL_SHA256 = ab3f0824e0d17bd9ce5075964350b0ab13bb8991e85a87f7d822be882b3a102e
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

fetch-git: .fetch-git-$(GIT_VERSION) .fetch-curl-$(CURL_VERSION)

.extract-git-$(GIT_VERSION): .fetch-git-$(GIT_VERSION)
	$(RM) $(GIT_BUILD_DIR)
	$(TAR) -xzf $(GIT_TARBALL)
	@$(TOUCH) $@

extract-git: .extract-git-$(GIT_VERSION)

.configure-git-$(GIT_VERSION): .extract-git-$(GIT_VERSION) .install-curl-$(CURL_VERSION)
	cd $(GIT_BUILD_DIR) && CFLAGS="-DCURL_STATICLIB" ./configure $(GIT_CONFIGURE_FLAGS) --prefix=$(GIT_INSTALL_DIR) --with-curl=$(CURL_INSTALL_DIR)
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

clean-git: clean-curl
	@$(RM) .build-git-$(GIT_VERSION)
	@$(RM) .configure-git-$(GIT_VERSION)
	@$(RM) .extract-git-$(GIT_VERSION)
	$(RM) $(GIT_BUILD_DIR)

distclean-git: clean-git distclean-curl
	@$(RM) .fetch-git-$(GIT_VERSION)
	$(RM) $(GIT_TARBALL)
	$(RM) $(GIT_MANPAGES_TARBALL)

env-git:
	@echo 'export PATH="$(GIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GIT_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## python-sphinx
##

PYTHON_SPHINX_VERSION = 1.2.3
PYTHON_SPHINX_TARBALL = Sphinx-$(PYTHON_SPHINX_VERSION).tar.gz
PYTHON_SPHINX_TARBALL_URL = http://pypi.python.org/packages/source/S/Sphinx/$(PYTHON_SPHINX_TARBALL)
PYTHON_SPHINX_TARBALL_SHA256 = 94933b64e2fe0807da0612c574a021c0dac28c7bd3c4a23723ae5a39ea8f3d04
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
	@echo 'export PATH="$(GRAPHVIZ_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GRAPHVIZ_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## GDB debugger
##

GDB_VERSION = 7.7.1
GDB_TARBALL_SHA256 = eefadb9831e3695d1eaef34e98b8f1fb441df6fe5071317ea49c6bd6ba213eff
GDB_TARBALL = gdb-$(GDB_VERSION).tar.gz
GDB_TARBALL_URL = http://ftp.gnu.org/gnu/gdb/$(GDB_TARBALL)
GDB_BUILD_DIR = gdb-$(GDB_VERSION)
GDB_INSTALL_DIR = $(PREFIX)/gdb-$(GDB_VERSION)

.fetch-gdb-$(GDB_VERSION):
	@$(RM) $(GDB_TARBALL)
	$(WGET) $(GDB_TARBALL_URL)
	@echo '$(GDB_TARBALL_SHA256)  $(GDB_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gdb: .fetch-gdb-$(GDB_VERSION)

.extract-gdb-$(GDB_VERSION): .fetch-gdb-$(GDB_VERSION)
	$(RM) $(GDB_BUILD_DIR)
	$(TAR) -xzf $(GDB_TARBALL)
	@$(TOUCH) $@

extract-gdb: .extract-gdb-$(GDB_VERSION)

.configure-gdb-$(GDB_VERSION): .extract-gdb-$(GDB_VERSION)
	cd $(GDB_BUILD_DIR) && ./configure $(GDB_CONFIGURE_FLAGS) --prefix=$(GDB_INSTALL_DIR)
	@$(TOUCH) $@

configure-gdb: .configure-gdb-$(GDB_VERSION)

.build-gdb-$(GDB_VERSION): .configure-gdb-$(GDB_VERSION)
	cd $(GDB_BUILD_DIR) && $(MAKE)
	@$(TOUCH) $@

build-gdb: .build-gdb-$(GDB_VERSION)

install-gdb: .build-gdb-$(GDB_VERSION)
	cd $(GDB_BUILD_DIR) && $(MAKE) install

clean-gdb:
	@$(RM) .build-gdb-$(GDB_VERSION)
	@$(RM) .extract-gdb-$(GDB_VERSION)
	@$(RM) .configure-gdb-$(GDB_VERSION)
	$(RM) $(GDB_BUILD_DIR)

distclean-gdb: clean-gdb
	@$(RM) .fetch-gdb-$(GDB_VERSION)
	$(RM) $(GDB_TARBALL)

env-gdb:
	@echo 'export PATH="$(GDB_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GDB_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(GDB_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Clang C++ compiler
##
CLANG_VERSION = 3.6.0
LLVM_TARBALL = llvm-$(CLANG_VERSION).src.tar.xz
LLVM_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(LLVM_TARBALL)
LLVM_TARBALL_SHA256 = b39a69e501b49e8f73ff75c9ad72313681ee58d6f430bfad4d81846fe92eb9ce
CFE_TARBALL = cfe-$(CLANG_VERSION).src.tar.xz
CFE_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(CFE_TARBALL)
CFE_TARBALL_SHA256 = be0e69378119fe26f0f2f74cffe82b7c26da840c9733fe522ed3c1b66b11082d
CLANG_BUILD_DIR = llvm-$(CLANG_VERSION).src
CLANG_CONFIGURE_FLAGS = --enable-optimized --enable-bindings=none --enable-shared
CLANG_BUILD_FLAGS = REQUIRES_RTTI=1
CLANG_INSTALL_DIR = $(PREFIX)/clang-$(CLANG_VERSION)

ifdef CLANG_GCC_TOOLCHAIN
    CLANG_CONFIGURE_FLAGS += --with-gcc-toolchain=$(CLANG_GCC_TOOLCHAIN)
endif

.fetch-clang-$(CLANG_VERSION):
	@$(RM) $(LLVM_TARBALL)
	@$(RM) $(CFE_TARBALL)
	$(WGET) $(LLVM_TARBALL_URL)
	$(WGET) $(CFE_TARBALL_URL)
	@echo '$(LLVM_TARBALL_SHA256)  $(LLVM_TARBALL)' | $(SHA256SUM)
	@echo '$(CFE_TARBALL_SHA256)  $(CFE_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-clang: .fetch-clang-$(CLANG_VERSION)

.extract-clang-$(CLANG_VERSION): .fetch-clang-$(CLANG_VERSION)
	$(RM) $(CLANG_BUILD_DIR)
	$(TAR) -xJf $(LLVM_TARBALL)
	cd $(CLANG_BUILD_DIR)/tools && $(TAR) -xJf $(CURDIR)/$(CFE_TARBALL) && mv cfe-$(CLANG_VERSION).src clang
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
	$(RM) $(CFE_TARBALL)

env-clang:
	@echo 'export PATH="$(CLANG_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(CLANG_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'
	@echo 'export LD_LIBRARY_PATH="$(CLANG_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'

##
## GMP (GNU Multiple Precision Arithmetic Library)
##

GMP_VERSION = 5.1.3
GMP_TARBALL = gmp-$(GMP_VERSION).tar.bz2
GMP_TARBALL_URL = http://ftp.gnu.org/gnu/gmp/$(GMP_TARBALL)
GMP_TARBALL_SHA256 = 752079520b4690531171d0f4532e40f08600215feefede70b24fabdc6f1ab160
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
	@$(RM) .install-gmp-$(GMP_VERSION)
	$(RM) $(GMP_TARBALL)
	$(RM) $(GMP_INSTALL_DIR)

##
## MPFR (Multiple-precision floating-point computations with correct rounding)
##

MPFR_VERSION = 3.1.2
MPFR_TARBALL = mpfr-$(MPFR_VERSION).tar.bz2
MPFR_TARBALL_URL = http://www.mpfr.org/mpfr-$(MPFR_VERSION)/$(MPFR_TARBALL)
MPFR_TARBALL_SHA256 = 79c73f60af010a30a5c27a955a1d2d01ba095b72537dab0ecaad57f5a7bb1b6b
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
	@$(RM) .install-mpfr-$(MPFR_VERSION)
	$(RM) $(MPFR_TARBALL)
	$(RM) $(MPFR_INSTALL_DIR)

##
## MPC (arithmetic of complex numbers with arbitrarily high precision and correct rounding)
##

MPC_VERSION = 1.0.3
MPC_TARBALL = mpc-$(MPC_VERSION).tar.gz
MPC_TARBALL_URL = http://www.multiprecision.org/mpc/download/$(MPC_TARBALL)
MPC_TARBALL_SHA256 = 617decc6ea09889fb08ede330917a00b16809b8db88c29c31bfbb49cbf88ecc3
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
	@$(RM) .install-mpc-$(MPC_VERSION)
	$(RM) $(MPC_TARBALL)
	$(RM) $(MPC_INSTALL_DIR)

##
## ISL (Integer Set Library)
##

ISL_VERSION = 0.14
ISL_TARBALL = isl-$(ISL_VERSION).tar.bz2
ISL_TARBALL_URL = ftp://ftp.fu-berlin.de/unix/languages/gcc/infrastructure/$(ISL_TARBALL)
ISL_TARBALL_SHA256 = 7e3c02ff52f8540f6a85534f54158968417fd676001651c8289c705bd0228f36
ISL_BUILD_DIR = isl-$(ISL_VERSION)
ISL_INSTALL_DIR = $(CURDIR)/.isl-$(ISL_VERSION)

.fetch-isl-$(ISL_VERSION):
	@$(RM) $(ISL_TARBALL)
	$(WGET) $(ISL_TARBALL_URL)
	@echo '$(ISL_TARBALL_SHA256)  $(ISL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-isl: .fetch-isl-$(ISL_VERSION)

.extract-isl-$(ISL_VERSION): .fetch-isl-$(ISL_VERSION)
	$(RM) $(ISL_BUILD_DIR)
	$(TAR) -xjf $(ISL_TARBALL)
	@$(TOUCH) $@

extract-isl: .extract-isl-$(ISL_VERSION)

.configure-isl-$(ISL_VERSION): .extract-isl-$(ISL_VERSION) .install-gmp-$(GMP_VERSION)
	cd $(ISL_BUILD_DIR) && LDFLAGS="-L$(GMP_INSTALL_DIR)/lib" ./configure --prefix=$(ISL_INSTALL_DIR) --libdir=$(ISL_INSTALL_DIR)/lib --with-gmp-prefix=$(GMP_INSTALL_DIR) --disable-shared
	@$(TOUCH) $@

configure-isl: .configure-isl-$(ISL_VERSION)

.build-isl-$(ISL_VERSION): .configure-isl-$(ISL_VERSION)
	cd $(ISL_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-isl: .build-isl-$(ISL_VERSION)

.install-isl-$(ISL_VERSION): .build-isl-$(ISL_VERSION)
	cd $(ISL_BUILD_DIR) && $(MAKE) install
	@$(TOUCH) $@

clean-isl:
	@$(RM) .build-isl-$(ISL_VERSION)
	@$(RM) .configure-isl-$(ISL_VERSION)
	@$(RM) .extract-isl-$(ISL_VERSION)
	$(RM) $(ISL_BUILD_DIR)

distclean-isl: clean-isl
	@$(RM) .fetch-isl-$(ISL_VERSION)
	@$(RM) .install-isl-$(ISL_VERSION)
	$(RM) $(ISL_TARBALL)
	$(RM) $(ISL_INSTALL_DIR)

##
## GCC (GNU Compiler Collection)
##

GCC_VERSION = 5.1.0
GCC_TARBALL = gcc-$(GCC_VERSION).tar.bz2
GCC_TARBALL_URL = http://ftp.gwdg.de/pub/misc/gcc/releases/gcc-$(GCC_VERSION)/$(GCC_TARBALL)
GCC_TARBALL_SHA256 = b7dafdf89cbb0e20333dbf5b5349319ae06e3d1a30bf3515b5488f7e89dca5ad
GCC_BUILD_DIR = gcc-$(GCC_VERSION)
GCC_BUILD_FLAGS = --enable-cxx-flags=-fPIC --enable-languages=c,c++,fortran,lto --disable-multilib
GCC_INSTALL_DIR = $(PREFIX)/gcc-$(GCC_VERSION)

.fetch-gcc-$(GCC_VERSION):
	@$(RM) $(GCC_TARBALL)
	$(WGET) $(GCC_TARBALL_URL)
	@echo '$(GCC_TARBALL_SHA256)  $(GCC_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gcc: .fetch-gcc-$(GCC_VERSION) .fetch-gmp-$(GMP_VERSION) .fetch-mpfr-$(MPFR_VERSION) .fetch-mpc-$(MPC_VERSION) .fetch-isl-$(ISL_VERSION)

.extract-gcc-$(GCC_VERSION): .fetch-gcc-$(GCC_VERSION) .extract-gmp-$(GMP_VERSION) .extract-mpfr-$(MPFR_VERSION) .extract-mpc-$(MPC_VERSION)
	$(RM) $(GCC_BUILD_DIR)
	$(TAR) -xjf $(GCC_TARBALL)
	$(CP) $(GMP_BUILD_DIR) $(GCC_BUILD_DIR)/gmp
	$(CP) $(MPFR_BUILD_DIR) $(GCC_BUILD_DIR)/mpfr
	$(CP) $(MPC_BUILD_DIR) $(GCC_BUILD_DIR)/mpc
	@$(TOUCH) $@

extract-gcc: .extract-gcc-$(GCC_VERSION)

.configure-gcc-$(GCC_VERSION): .extract-gcc-$(GCC_VERSION) .install-isl-$(ISL_VERSION)
	cd $(GCC_BUILD_DIR) && ac_cv_lib_pwl_PWL_handle_timeout=no ./configure --prefix=$(GCC_INSTALL_DIR) --with-isl=$(ISL_INSTALL_DIR) --enable-build-with-cxx $(GCC_BUILD_FLAGS)
	@$(TOUCH) $@

configure-gcc: .configure-gcc-$(GCC_VERSION)

.build-gcc-$(GCC_VERSION): .configure-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gcc: .build-gcc-$(GCC_VERSION)

install-gcc: .build-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && $(MAKE) install
	cd $(GMP_INSTALL_DIR) && $(CP) include `$(GCC_INSTALL_DIR)/bin/gcc -print-file-name=plugin`

clean-gcc: clean-gmp clean-mpfr clean-mpc clean-isl
	@$(RM) .build-gcc-$(GCC_VERSION)
	@$(RM) .configure-gcc-$(GCC_VERSION)
	@$(RM) .extract-gcc-$(GCC_VERSION)
	$(RM) $(GCC_BUILD_DIR)

distclean-gcc: clean-gcc distclean-gmp distclean-mpfr distclean-mpc distclean-isl
	@$(RM) .fetch-gcc-$(GCC_VERSION)
	$(RM) $(GCC_TARBALL)

env-gcc:
	@echo 'export PATH="$(GCC_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GCC_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'
	@echo 'export LD_LIBRARY_PATH="$(GCC_INSTALL_DIR)/lib64$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'

##
## HALMD Highly Accerelated Large-scale Molecular Dynamics
##

HALMD_VERSION = 1.0-alpha4
HALMD_TARBALL = halmd-$(HALMD_VERSION).tar.bz2
HALMD_TARBALL_URL = http://code.halmd.org/tar/$(HALMD_TARBALL)
HALMD_SOURCE_DIR = halmd-$(HALMD_VERSION)
HALMD_BUILD_DIR = $(HALMD_SOURCE_DIR)/build/release
HALMD_INSTALL_DIR = $(PREFIX)/halmd-$(HALMD_VERSION)
HALMD_BUILD_ENV =

define HALMD_CERT_PUBKEY
LS0tLS1CRUdJTiBQVUJMSUMgS0VZLS0tLS0KTUlJQklqQU5CZ2txaGtpRzl3MEJBUUVGQUFPQ0FR
OEFNSUlCQ2dLQ0FRRUEweFoycHRJYjNYOVdaYTNQWDhmVgpsVUtBaGpHUVA0YUYrZzk0Mkh0OXp0
SDZFM2NhN3lJUzdYdk9Dd2RpVFpuQy9WN05KNEw4NloxRkI4bmYrajEwCjlyaEJBWWkyaVZIZzNo
WXBiMldtWjRCamJ3YWpDdGRPS2pYNzZxUFlYay8zTDVGSlZFSUM3S1hUOFpGMmVVTnkKcFlsZWN5
T1RFSjREVUJZS2VsZEFzZlgwRGMvSCtabjdmK1lBem54cStZTzZQdXUzU3dweVhFNkd5eGNHWlNW
cgpEZ2ZYM3o4bEM1WG1aUjZ6WUFYOU9qUzFubnVjNGJCQUt4S3RUWFh2eVZvMTNodU0rbU1kUlZJ
WUlnMFMzMGV1CnFxQW0zbnc0MlRZM0d2c1dzU0h0RVA3SGR0Q1ZlaUJvQ2l4SlN1Q0lGR1c3T2Vm
MzJEZ1VndnBObHl3RlkrdHUKeFFJREFRQUIKLS0tLS1FTkQgUFVCTElDIEtFWS0tLS0tCg==
endef
export HALMD_CERT_PUBKEY

.fetch-halmd-$(HALMD_VERSION):
	$(RM) $(HALMD_TARBALL)
	$(WGET) $(HALMD_TARBALL_URL)
	$(WGET) -q $(HALMD_TARBALL_URL).sig
	echo $$HALMD_CERT_PUBKEY | $(BASE64) | $(OPENSSL) dgst -sha512 -verify /dev/stdin \
	  -signature $(HALMD_TARBALL).sig $(HALMD_TARBALL)
	$(RM) $(HALMD_TARBALL).sig
	@$(TOUCH) $@

fetch-halmd: .fetch-halmd-$(HALMD_VERSION)

.extract-halmd-$(HALMD_VERSION): .fetch-halmd-$(HALMD_VERSION)
	$(TAR) -xjf $(HALMD_TARBALL)
	@$(TOUCH) $@

extract-halmd: .extract-halmd-$(HALMD_VERSION)

.configure-halmd-$(HALMD_VERSION): .extract-halmd-$(HALMD_VERSION)
	$(MKDIR) $(HALMD_BUILD_DIR)
	cd $(HALMD_BUILD_DIR) && $(HALMD_BUILD_ENV) $(CMAKE) -DCMAKE_INSTALL_PREFIX=$(HALMD_INSTALL_DIR) $(CURDIR)/$(HALMD_SOURCE_DIR)
	@$(TOUCH) $@

configure-halmd: .configure-halmd-$(HALMD_VERSION)

.build-halmd-$(HALMD_VERSION): .configure-halmd-$(HALMD_VERSION)
	cd $(HALMD_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-halmd: .build-halmd-$(HALMD_VERSION)

test-halmd: .build-halmd-$(HALMD_VERSION)
	cd $(HALMD_BUILD_DIR) && $(CTEST)

install-halmd: .build-halmd-$(HALMD_VERSION)
	cd $(HALMD_BUILD_DIR) && $(MAKE) install

clean-halmd:
	@$(RM) .build-halmd-$(HALMD_VERSION)
	@$(RM) .configure-halmd-$(HALMD_VERSION)
	$(RM) $(HALMD_BUILD_DIR)

distclean-halmd: clean-halmd
	@$(RM) .fetch-halmd-$(HALMD_VERSION)
	@$(RM) .extract-halmd-$(HALMD_VERSION)
	$(RM) $(HALMD_SOURCE_DIR)
	$(RM) $(HALMD_TARBALL)

env-halmd:
	@echo
	@echo 'export PATH="$(HALMD_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(HALMD_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## nvCUDA-tools - A collection of tools for NVIDIA CUDA compute devices
##

NVCUDA_TOOLS_VERSION = master
NVCUDA_TOOLS_GIT_URL = http://github.com/halmd-org/nvcuda-tools.git
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
	cd $(NVCUDA_TOOLS_BUILD_DIR) && $(CMAKE) -DCMAKE_INSTALL_PREFIX=$(NVCUDA_TOOLS_INSTALL_DIR) $(CURDIR)/$(NVCUDA_TOOLS_SOURCE_DIR)
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
	@echo 'export PATH="$(NVCUDA_TOOLS_INSTALL_DIR)/bin$${PATH+:$$PATH}"'

##
## Ninja - a small build system with a focus on speed
##
NINJA_VERSION = 1.5.3
NINJA_TARBALL = ninja-$(NINJA_VERSION).tar.gz
NINJA_TARBALL_REMOTE = v$(NINJA_VERSION).tar.gz
NINJA_TARBALL_URL = https://github.com/martine/ninja/archive/$(NINJA_TARBALL_REMOTE)
NINJA_TARBALL_SHA256 = 7c953b5a7c26cfcd082882e3f3e2cd08fee8848ad228bb47223b18ea18777ec0
NINJA_BUILD_DIR = ninja-$(NINJA_VERSION)
NINJA_CONFIGURE_FLAGS =
NINJA_INSTALL_DIR = $(PREFIX)/ninja-$(NINJA_VERSION)

.fetch-ninja-$(NINJA_VERSION):
	@$(RM) $(NINJA_TARBALL)
	$(WGET) $(NINJA_TARBALL_URL) -O $(NINJA_TARBALL)
	@echo '$(NINJA_TARBALL_SHA256)  $(NINJA_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-ninja: .fetch-ninja-$(NINJA_VERSION)

.extract-ninja-$(NINJA_VERSION): .fetch-ninja-$(NINJA_VERSION)
	$(RM) $(NINJA_BUILD_DIR)
	$(TAR) -xzf $(NINJA_TARBALL)
	@$(TOUCH) $@

extract-ninja: .extract-ninja-$(NINJA_VERSION)

.build-ninja-$(NINJA_VERSION): .extract-ninja-$(NINJA_VERSION)
	cd $(NINJA_BUILD_DIR) && ./configure.py --bootstrap $(NINJA_CONFIGURE_FLAGS)
	@$(TOUCH) $@

build-ninja: .build-ninja-$(NINJA_VERSION)

install-ninja: .build-ninja-$(NINJA_VERSION)
	install -d $(NINJA_INSTALL_DIR)/bin
	cd $(NINJA_BUILD_DIR) && cp ninja $(NINJA_INSTALL_DIR)/bin

clean-ninja:
	@$(RM) .build-ninja-$(NINJA_VERSION)
	@$(RM) .extract-ninja-$(NINJA_VERSION)
	$(RM) $(NINJA_BUILD_DIR)

distclean-ninja: clean-ninja
	@$(RM) .fetch-ninja-$(NINJA_VERSION)
	$(RM) $(NINJA_TARBALL)

env-ninja:
	@echo
	@echo 'export PATH="$(NINJA_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
