#
# Copyright © 2011-2017 Peter Colberg
# Copyright © 2013-2025 Felix Höfling
# Copyright © 2013      Nicolas Höft
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
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
LN        = ln -sf
OPENSSL   = openssl
PATCH     = patch
RM        = rm -rf
SED       = sed
SHA256SUM = sha256sum --check
SHA512SUM = sha512sum --check
TAR       = tar
TOUCH     = touch
WGET      = wget

##
## define top-level targets
##

build: build-cmake build-lua build-luajit build-boost build-hdf5 build-git build-python-sphinx build-graphviz build-gdb build-llvm build-clang build-gcc build-halmd build-nvcuda-tools build-ninja

fetch: fetch-cmake fetch-lua fetch-luajit fetch-boost fetch-hdf5 fetch-git fetch-python-sphinx fetch-graphviz fetch-gdb fetch-llvm fetch-clang fetch-gcc fetch-halmd fetch-nvcuda-tools fetch-ninja

install: install-cmake install-lua install-luajit install-boost install-hdf5 install-git install-python-sphinx install-graphviz install-gdb install-llvm install-clang install-gcc install-halmd install-nvcuda-tools install-ninja

clean: clean-cmake clean-lua clean-luajit clean-boost clean-hdf5 clean-git clean-python-sphinx clean-graphviz clean-gdb clean-llvm clean-clang clean-gcc clean-halmd clean-nvcuda-tools clean-ninja

distclean: distclean-cmake distclean-lua distclean-luajit distclean-boost distclean-hdf5 distclean-git distclean-python-sphinx distclean-graphviz distclean-gdb distclean-llvm distclean-clang distclean-gcc distclean-halmd distclean-nvcuda-tools distclean-ninja

env: env-cmake env-lua env-luajit env-boost env-hdf5 env-git env-python-sphinx env-graphviz env-gdb env-llvm env-clang env-gcc env-halmd env-nvcuda-tools env-ninja

##
## CMake
##

CMAKE_MAJOR_VERSION = 3
CMAKE_MINOR_VERSION = 31
CMAKE_PATCH_VERSION = 4
CMAKE_VERSION = $(CMAKE_MAJOR_VERSION).$(CMAKE_MINOR_VERSION).$(CMAKE_PATCH_VERSION)
CMAKE_TARBALL = cmake-$(CMAKE_VERSION).tar.gz
CMAKE_TARBALL_URL = https://cmake.org/files/v$(CMAKE_MAJOR_VERSION).$(CMAKE_MINOR_VERSION)/$(CMAKE_TARBALL)
CMAKE_TARBALL_SHA256 = a6130bfe75f5ba5c73e672e34359f7c0a1931521957e8393a5c2922c8b0f7f25
CMAKE_SOURCE_DIR = cmake-$(CMAKE_VERSION)
CMAKE_BUILD_DIR = $(CMAKE_SOURCE_DIR)/build
CMAKE_CONFIGURE_FLAGS = --sphinx-man
CMAKE_INSTALL_DIR = $(PREFIX)/cmake-$(CMAKE_VERSION)

.fetch-cmake-$(CMAKE_VERSION):
	@$(RM) $(CMAKE_TARBALL)
	$(WGET) $(CMAKE_TARBALL_URL) -O $(CMAKE_TARBALL)
	@echo '$(CMAKE_TARBALL_SHA256)  $(CMAKE_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-cmake: .fetch-cmake-$(CMAKE_VERSION)

.extract-cmake-$(CMAKE_VERSION): .fetch-cmake-$(CMAKE_VERSION)
	$(RM) $(CMAKE_SOURCE_DIR)
	$(TAR) -xzf $(CMAKE_TARBALL)
	@$(TOUCH) $@

extract-cmake: .extract-cmake-$(CMAKE_VERSION)

.configure-cmake-$(CMAKE_VERSION): .extract-cmake-$(CMAKE_VERSION)
	$(MKDIR) $(CMAKE_BUILD_DIR)
	cd $(CMAKE_BUILD_DIR) && ../configure --prefix=$(CMAKE_INSTALL_DIR) $(CMAKE_CONFIGURE_FLAGS)
	@$(TOUCH) $@

configure-cmake: .configure-cmake-$(CMAKE_VERSION)

.build-cmake-$(CMAKE_VERSION): .configure-cmake-$(CMAKE_VERSION)
	cd $(CMAKE_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-cmake: .build-cmake-$(CMAKE_VERSION)

install-cmake: .build-cmake-$(CMAKE_VERSION)
	cd $(CMAKE_BUILD_DIR) && $(MAKE) install

clean-cmake:
	@$(RM) .build-cmake-$(CMAKE_VERSION)
	@$(RM) .configure-cmake-$(CMAKE_VERSION)
	$(RM) $(CMAKE_BUILD_DIR)

distclean-cmake: clean-cmake
	@$(RM) .fetch-cmake-$(CMAKE_VERSION)
	@$(RM) .extract-cmake-$(CMAKE_VERSION)
	$(RM) $(CMAKE_SOURCE_DIR)
	$(RM) $(CMAKE_TARBALL)

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
LUA_VERSION = 5.2.4
LUA_TARBALL_SHA256 = b9e2e4aad6789b3b63a056d442f7b39f0ecfca3ae0f1fc0ae4e9614401b69f4b
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

LUAJIT_VERSION = 2.1
LUAJIT_ARCHIVE = LuaJIT-$(LUAJIT_VERSION).tar.gz
LUAJIT_ARCHIVE_URL = https://github.com/LuaJIT/LuaJIT/archive/refs/heads/v$(LUAJIT_VERSION).tar.gz
LUAJIT_BUILD_DIR = LuaJIT-$(LUAJIT_VERSION)
LUAJIT_INSTALL_DIR = $(PREFIX)/luajit-$(LUAJIT_VERSION)
LUAJIT_CFLAGS = -fPIC -DLUAJIT_ENABLE_LUA52COMPAT -DLUAJIT_CPU_SSE2

ifdef USE_VALGRIND
LUAJIT_CFLAGS += -DLUAJIT_USE_VALGRIND
endif

.fetch-luajit-$(LUAJIT_VERSION):
	@$(RM) $(LUAJIT_ARCHIVE)
	$(WGET) $(LUAJIT_ARCHIVE_URL) -O $(LUAJIT_ARCHIVE)
	@$(TOUCH) $@

fetch-luajit: .fetch-luajit-$(LUAJIT_VERSION)

.extract-luajit-$(LUAJIT_VERSION): .fetch-luajit-$(LUAJIT_VERSION)
	$(RM) $(LUAJIT_BUILD_DIR)
	$(TAR) -xzf $(LUAJIT_ARCHIVE)
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
	$(RM) $(LUAJIT_ARCHIVE)

env-luajit:
	@echo 'export PATH="$(LUAJIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(LUAJIT_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LUAJIT_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## luatrace
##

LUATRACE_VERSION = d9d8918
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

BOOST_VERSION = 1.87.0
BOOST_RELEASE = 1_87_0
BOOST_ABI = c++11
BOOST_TOOLSET = gcc
BOOST_TARBALL = boost_$(BOOST_RELEASE).tar.bz2
BOOST_TARBALL_URL = https://archives.boost.io/release/$(BOOST_VERSION)/source/$(BOOST_TARBALL)
BOOST_TARBALL_SHA256 = af57be25cb4c4f4b413ed692fe378affb4352ea50fbe294a11ef548f4d527d89
BOOST_BUILD_DIR = boost_$(BOOST_RELEASE)
BOOST_INSTALL_DIR = $(PREFIX)/boost_$(BOOST_RELEASE)
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
	cd $(BOOST_BUILD_DIR) && ./b2 $(BOOST_BUILD_FLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-boost: .build-boost-$(BOOST_VERSION)

install-boost: .build-boost-$(BOOST_VERSION)
	cd $(BOOST_BUILD_DIR) && ./b2 $(BOOST_BUILD_FLAGS) install --prefix=$(BOOST_INSTALL_DIR)

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
HDF5_VERSION = 1.14.5
HDF5_MAJOR_VERSION = $(shell echo $(HDF5_VERSION) | cut -f1 -d.)
HDF5_MINOR_VERSION = $(shell echo $(HDF5_VERSION) | cut -f2 -d.)
HDF5_TARBALL = hdf5-$(HDF5_VERSION).tar.bz2
HDF5_TARBALL_URL = https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_MAJOR_VERSION).$(HDF5_MINOR_VERSION)/hdf5-$(HDF5_VERSION)/src/$(HDF5_TARBALL)
HDF5_TARBALL_SHA256 = ec2e13c52e60f9a01491bb3158cb3778c985697131fc6a342262d32a26e58e44
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
GIT_VERSION = 2.18.0
GIT_TARBALL = git-$(GIT_VERSION).tar.gz
GIT_TARBALL_URL = https://www.kernel.org/pub/software/scm/git/$(GIT_TARBALL)
GIT_TARBALL_SHA256 = 94faf2c0b02a7920b0b46f4961d8e9cad08e81418614102898a55f980fa3e7e4
GIT_MANPAGES_TARBALL = git-manpages-$(GIT_VERSION).tar.gz
GIT_MANPAGES_TARBALL_URL = https://www.kernel.org/pub/software/scm/git/$(GIT_MANPAGES_TARBALL)
GIT_MANPAGES_TARBALL_SHA256 = 6cf38ab3ad43ccdcd6a73ffdcf2a016d56ab6b4b240a574b0bb96f520a04ff55
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

PYTHON_SPHINX_VERSION = 8.1.3
PYTHON_SPHINX_TARBALL = Sphinx-$(PYTHON_SPHINX_VERSION).tar.gz
PYTHON_SPHINX_TARBALL_URL = https://github.com/sphinx-doc/sphinx/archive/v${PYTHON_SPHINX_VERSION}.tar.gz
PYTHON_SPHINX_TARBALL_SHA256 = 0fcc28999fe8e4fcc49a4ab01e3e987f6fbb3af32995db74e6fc8f8d01dcaaca
PYTHON_SPHINX_BUILD_DIR = Sphinx-$(PYTHON_SPHINX_VERSION)
PYTHON_SPHINX_INSTALL_DIR = $(PREFIX)/python-sphinx-$(PYTHON_SPHINX_VERSION)
PYTHON_SPHINX_PYTHONPATH = $(PYTHON_SPHINX_INSTALL_DIR)/lib/python

.fetch-python-sphinx-$(PYTHON_SPHINX_VERSION):
	@$(RM) $(PYTHON_SPHINX_TARBALL)
	$(WGET) $(PYTHON_SPHINX_TARBALL_URL) -O $(PYTHON_SPHINX_TARBALL)
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
## LLVM compiler infrastructure
##

LLVM_VERSION = 13.0.1
LLVM_TARBALL = llvm-$(LLVM_VERSION).src.tar.xz
LLVM_TARBALL_URL = https://github.com/llvm/llvm-project/releases/download/llvmorg-$(LLVM_VERSION)/$(LLVM_TARBALL)
LLVM_TARBALL_SHA256 = ec6b80d82c384acad2dc192903a6cf2cdbaffb889b84bfb98da9d71e630fc834
LLVM_SOURCE_DIR = llvm-$(LLVM_VERSION).src
LLVM_BUILD_DIR = $(LLVM_SOURCE_DIR)/build
LLVM_CONFIGURE_FLAGS = -DCMAKE_BUILD_TYPE=Release -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_ENABLE_RTTI=ON -DLLVM_BUILD_LLVM_DYLIB=ON -DLLVM_LINK_LLVM_DYLIB=ON
LLVM_INSTALL_DIR = $(PREFIX)/llvm-$(LLVM_VERSION)

.fetch-llvm-$(LLVM_VERSION):
	@$(RM) $(LLVM_TARBALL)
	$(WGET) -O $(LLVM_TARBALL) $(LLVM_TARBALL_URL)
	@echo '$(LLVM_TARBALL_SHA256)  $(LLVM_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-llvm: .fetch-llvm-$(LLVM_VERSION)

.extract-llvm-$(LLVM_VERSION): .fetch-llvm-$(LLVM_VERSION)
	$(RM) $(LLVM_SOURCE_DIR)
	$(TAR) -xf $(LLVM_TARBALL)
	@$(TOUCH) $@

extract-llvm: .extract-llvm-$(LLVM_VERSION)

.configure-llvm-$(LLVM_VERSION): .extract-llvm-$(LLVM_VERSION)
	$(MKDIR) $(LLVM_BUILD_DIR)
	cd $(LLVM_BUILD_DIR) && $(CMAKE) -DCMAKE_INSTALL_PREFIX=$(LLVM_INSTALL_DIR) $(LLVM_CONFIGURE_FLAGS) ..
	@$(TOUCH) $@

configure-llvm: .configure-llvm-$(LLVM_VERSION)

.build-llvm-$(LLVM_VERSION): .configure-llvm-$(LLVM_VERSION)
	cd $(LLVM_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-llvm: .build-llvm-$(LLVM_VERSION)

install-llvm: .build-llvm-$(LLVM_VERSION)
	cd $(LLVM_BUILD_DIR) && $(MAKE) install

clean-llvm:
	@$(RM) .build-llvm-$(LLVM_VERSION)
	@$(RM) .configure-llvm-$(LLVM_VERSION)
	@$(RM) .extract-llvm-$(LLVM_VERSION)
	$(RM) $(LLVM_SOURCE_DIR)

distclean-llvm: clean-llvm
	@$(RM) .fetch-llvm-$(LLVM_VERSION)
	$(RM) $(LLVM_TARBALL)

env-llvm:
	@echo 'export PATH="$(LLVM_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export CPATH="$(LLVM_INSTALL_DIR)/include$${CPATH+:$$CPATH}"'
	@echo 'export LIBRARY_PATH="$(LLVM_INSTALL_DIR)/lib$${LIBRARY_PATH+:$$LIBRARY_PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(LLVM_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LLVM_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Clang C++ compiler
##

CLANG_VERSION = 13.0.1
CLANG_TARBALL = clang-$(CLANG_VERSION).src.tar.xz
CLANG_TARBALL_URL = https://github.com/llvm/llvm-project/releases/download/llvmorg-$(CLANG_VERSION)/$(CLANG_TARBALL)
CLANG_TARBALL_SHA256 = 787a9e2d99f5c8720aa1773e4be009461cd30d3bd40fdd24591e473467c917c9
CLANG_SOURCE_DIR = clang-$(CLANG_VERSION).src
CLANG_BUILD_DIR = $(CLANG_SOURCE_DIR)/build
CLANG_CONFIGURE_FLAGS = -DCMAKE_BUILD_TYPE=Release
CLANG_INSTALL_DIR = $(PREFIX)/clang-$(CLANG_VERSION)

ifndef USE_SYSTEM_GCC
CLANG_CONFIGURE_FLAGS += -DGCC_INSTALL_PREFIX=$(GCC_INSTALL_DIR)
endif

.fetch-clang-$(CLANG_VERSION):
	@$(RM) $(CLANG_TARBALL)
	$(WGET) -O $(CLANG_TARBALL) $(CLANG_TARBALL_URL)
	@echo '$(CLANG_TARBALL_SHA256)  $(CLANG_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-clang: .fetch-clang-$(CLANG_VERSION)

.extract-clang-$(CLANG_VERSION): .fetch-clang-$(CLANG_VERSION)
	$(RM) $(CLANG_SOURCE_DIR)
	$(TAR) -xf $(CLANG_TARBALL)
	@$(TOUCH) $@

extract-clang: .extract-clang-$(CLANG_VERSION)

.configure-clang-$(CLANG_VERSION): .extract-clang-$(CLANG_VERSION)
	$(MKDIR) $(CLANG_BUILD_DIR)
	cd $(CLANG_BUILD_DIR) && $(CMAKE) -DCMAKE_INSTALL_PREFIX=$(CLANG_INSTALL_DIR) $(CLANG_CONFIGURE_FLAGS) ..
	@$(TOUCH) $@

configure-clang: .configure-clang-$(CLANG_VERSION)

.build-clang-$(CLANG_VERSION): .configure-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-clang: .build-clang-$(CLANG_VERSION)

install-clang: .build-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && $(MAKE) install

clean-clang:
	@$(RM) .build-clang-$(CLANG_VERSION)
	@$(RM) .configure-clang-$(CLANG_VERSION)
	@$(RM) .extract-clang-$(CLANG_VERSION)
	$(RM) $(CLANG_SOURCE_DIR)

distclean-clang: clean-clang
	@$(RM) .fetch-clang-$(CLANG_VERSION)
	$(RM) $(CLANG_TARBALL)

env-clang:
	@echo 'export PATH="$(CLANG_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(CLANG_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## GMP (GNU Multiple Precision Arithmetic Library)
##

GMP_VERSION = 6.2.1
GMP_TARBALL = gmp-$(GMP_VERSION).tar.bz2
GMP_TARBALL_URL = ftp://ftp.fu-berlin.de/unix/languages/gcc/infrastructure/$(GMP_TARBALL)
GMP_TARBALL_SHA512 = 8904334a3bcc5c896ececabc75cda9dec642e401fb5397c4992c4fabea5e962c9ce8bd44e8e4233c34e55c8010cc28db0545f5f750cbdbb5f00af538dc763be9
GMP_BUILD_DIR = gmp-$(GMP_VERSION)

.fetch-gmp-$(GMP_VERSION):
	@$(RM) $(GMP_TARBALL)
	$(WGET) $(GMP_TARBALL_URL)
	@echo '$(GMP_TARBALL_SHA512)  $(GMP_TARBALL)' | $(SHA512SUM)
	@$(TOUCH) $@

fetch-gmp: .fetch-gmp-$(GMP_VERSION)

.extract-gmp-$(GMP_VERSION): .fetch-gmp-$(GMP_VERSION)
	$(RM) $(GMP_BUILD_DIR)
	$(TAR) -xjf $(GMP_TARBALL)
	@$(TOUCH) $@

extract-gmp: .extract-gmp-$(GMP_VERSION)

clean-gmp:
	@$(RM) .extract-gmp-$(GMP_VERSION)
	$(RM) $(GMP_BUILD_DIR)

distclean-gmp: clean-gmp
	@$(RM) .fetch-gmp-$(GMP_VERSION)
	$(RM) $(GMP_TARBALL)

##
## MPFR (Multiple-precision floating-point computations with correct rounding)
##

MPFR_VERSION = 4.1.0
MPFR_TARBALL = mpfr-$(MPFR_VERSION).tar.bz2
MPFR_TARBALL_URL = ftp://ftp.fu-berlin.de/unix/languages/gcc/infrastructure/$(MPFR_TARBALL)
MPFR_TARBALL_SHA512 = 410208ee0d48474c1c10d3d4a59decd2dfa187064183b09358ec4c4666e34d74383128436b404123b831e585d81a9176b24c7ced9d913967c5fce35d4040a0b4
MPFR_BUILD_DIR = mpfr-$(MPFR_VERSION)

.fetch-mpfr-$(MPFR_VERSION):
	@$(RM) $(MPFR_TARBALL)
	$(WGET) $(MPFR_TARBALL_URL)
	@echo '$(MPFR_TARBALL_SHA512)  $(MPFR_TARBALL)' | $(SHA512SUM)
	@$(TOUCH) $@

fetch-mpfr: .fetch-mpfr-$(MPFR_VERSION)

.extract-mpfr-$(MPFR_VERSION): .fetch-mpfr-$(MPFR_VERSION)
	$(RM) $(MPFR_BUILD_DIR)
	$(TAR) -xjf $(MPFR_TARBALL)
	@$(TOUCH) $@

extract-mpfr: .extract-mpfr-$(MPFR_VERSION)

clean-mpfr:
	@$(RM) .extract-mpfr-$(MPFR_VERSION)
	$(RM) $(MPFR_BUILD_DIR)

distclean-mpfr: clean-mpfr
	@$(RM) .fetch-mpfr-$(MPFR_VERSION)
	$(RM) $(MPFR_TARBALL)

##
## MPC (arithmetic of complex numbers with arbitrarily high precision and correct rounding)
##

MPC_VERSION = 1.2.1
MPC_TARBALL = mpc-$(MPC_VERSION).tar.gz
MPC_TARBALL_URL = ftp://ftp.fu-berlin.de/unix/languages/gcc/infrastructure/$(MPC_TARBALL)
MPC_TARBALL_SHA512 = 3279f813ab37f47fdcc800e4ac5f306417d07f539593ca715876e43e04896e1d5bceccfb288ef2908a3f24b760747d0dbd0392a24b9b341bc3e12082e5c836ee
MPC_BUILD_DIR = mpc-$(MPC_VERSION)

.fetch-mpc-$(MPC_VERSION):
	@$(RM) $(MPC_TARBALL)
	$(WGET) $(MPC_TARBALL_URL)
	@echo '$(MPC_TARBALL_SHA512)  $(MPC_TARBALL)' | $(SHA512SUM)
	@$(TOUCH) $@

fetch-mpc: .fetch-mpc-$(MPC_VERSION)

.extract-mpc-$(MPC_VERSION): .fetch-mpc-$(MPC_VERSION)
	$(RM) $(MPC_BUILD_DIR)
	$(TAR) -xzf $(MPC_TARBALL)
	@$(TOUCH) $@

extract-mpc: .extract-mpc-$(MPC_VERSION)

clean-mpc:
	@$(RM) .extract-mpc-$(MPC_VERSION)
	$(RM) $(MPC_BUILD_DIR)

distclean-mpc: clean-mpc
	@$(RM) .fetch-mpc-$(MPC_VERSION)
	$(RM) $(MPC_TARBALL)

##
## ISL (Integer Set Library)
##

ISL_VERSION = 0.24
ISL_TARBALL = isl-$(ISL_VERSION).tar.bz2
ISL_TARBALL_URL = ftp://ftp.fu-berlin.de/unix/languages/gcc/infrastructure/$(ISL_TARBALL)
ISL_TARBALL_SHA512 = aab3bddbda96b801d0f56d2869f943157aad52a6f6e6a61745edd740234c635c38231af20bc3f1a08d416a5e973a90e18249078ed8e4ae2f1d5de57658738e95
ISL_BUILD_DIR = isl-$(ISL_VERSION)

.fetch-isl-$(ISL_VERSION):
	@$(RM) $(ISL_TARBALL)
	$(WGET) $(ISL_TARBALL_URL)
	@echo '$(ISL_TARBALL_SHA512)  $(ISL_TARBALL)' | $(SHA512SUM)
	@$(TOUCH) $@

fetch-isl: .fetch-isl-$(ISL_VERSION)

.extract-isl-$(ISL_VERSION): .fetch-isl-$(ISL_VERSION)
	$(RM) $(ISL_BUILD_DIR)
	$(TAR) -xjf $(ISL_TARBALL)
	@$(TOUCH) $@

extract-isl: .extract-isl-$(ISL_VERSION)

clean-isl:
	@$(RM) .extract-isl-$(ISL_VERSION)
	$(RM) $(ISL_BUILD_DIR)

distclean-isl: clean-isl
	@$(RM) .fetch-isl-$(ISL_VERSION)
	$(RM) $(ISL_TARBALL)

##
## GCC (GNU Compiler Collection)
##

GCC_VERSION = 14.2.0
GCC_TARBALL = gcc-$(GCC_VERSION).tar.xz
GCC_TARBALL_URL = http://ftp.gwdg.de/pub/misc/gcc/releases/gcc-$(GCC_VERSION)/$(GCC_TARBALL)
GCC_TARBALL_SHA512 = d6828a5702ff4b667cc3e1e7e9f180191041b7efb68ffdc54248a42aa1799f41db6743acfe9ab74ea59977ba06f425fcf943a9fe3a77f9db706fc6bdbd657c1a
GCC_BUILD_DIR = gcc-$(GCC_VERSION)
GCC_BUILD_FLAGS = --enable-cxx-flags=-fPIC --enable-languages=c,c++,fortran,lto --disable-multilib
GCC_INSTALL_DIR = $(PREFIX)/gcc-$(GCC_VERSION)

.fetch-gcc-$(GCC_VERSION):
	@$(RM) $(GCC_TARBALL)
	$(WGET) $(GCC_TARBALL_URL)
	@echo '$(GCC_TARBALL_SHA512)  $(GCC_TARBALL)' | $(SHA512SUM)
	@$(TOUCH) $@

fetch-gcc: .fetch-gcc-$(GCC_VERSION) .fetch-gmp-$(GMP_VERSION) .fetch-mpfr-$(MPFR_VERSION) .fetch-mpc-$(MPC_VERSION) .fetch-isl-$(ISL_VERSION)

.extract-gcc-$(GCC_VERSION): .fetch-gcc-$(GCC_VERSION) .extract-gmp-$(GMP_VERSION) .extract-mpfr-$(MPFR_VERSION) .extract-mpc-$(MPC_VERSION) .extract-isl-$(ISL_VERSION)
	$(RM) $(GCC_BUILD_DIR)
	$(TAR) -xJf $(GCC_TARBALL)
	$(LN) ../$(GMP_BUILD_DIR) $(GCC_BUILD_DIR)/gmp
	$(LN) ../$(MPFR_BUILD_DIR) $(GCC_BUILD_DIR)/mpfr
	$(LN) ../$(MPC_BUILD_DIR) $(GCC_BUILD_DIR)/mpc
	$(LN) ../$(ISL_BUILD_DIR) $(GCC_BUILD_DIR)/isl
	@$(TOUCH) $@

extract-gcc: .extract-gcc-$(GCC_VERSION)

.configure-gcc-$(GCC_VERSION): .extract-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && ac_cv_lib_pwl_PWL_handle_timeout=no ./configure --prefix=$(GCC_INSTALL_DIR) --enable-build-with-cxx $(GCC_BUILD_FLAGS)
	@$(TOUCH) $@

configure-gcc: .configure-gcc-$(GCC_VERSION)

.build-gcc-$(GCC_VERSION): .configure-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && $(MAKE) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gcc: .build-gcc-$(GCC_VERSION)

install-gcc: .build-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && $(MAKE) install

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

HALMD_VERSION = 1.0.0
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
NINJA_VERSION = 1.12.1
NINJA_TARBALL = ninja-$(NINJA_VERSION).tar.gz
NINJA_TARBALL_REMOTE = v$(NINJA_VERSION).tar.gz
NINJA_TARBALL_URL = https://github.com/ninja-build/ninja/archive/$(NINJA_TARBALL_REMOTE)
NINJA_TARBALL_SHA256 = 821bdff48a3f683bc4bb3b6f0b5fe7b2d647cf65d52aeb63328c91a6c6df285a
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
