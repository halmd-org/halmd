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

WGET = wget
TAR = tar
UNZIP = unzip
RM = rm -rf
CP = cp -r
TOUCH = touch
PATCH = patch
SHA256SUM = sha256sum --check

##
## define top-level targets
##

build: build-cmake build-lua build-boost build-luabind build-hdf5

fetch: fetch-cmake fetch-lua fetch-boost fetch-luabind fetch-hdf5

install: install-cmake install-lua install-boost install-luabind install-hdf5

clean: clean-cmake clean-lua clean-boost clean-luabind clean-hdf5

distclean: distclean-cmake distclean-lua distclean-boost distclean-luabind distclean-hdf5

env: env-cmake env-lua env-boost env-luabind env-hdf5

##
## CMake with CMake-CUDA patch
##

CMAKE_VERSION = 2.8.8
CMAKE_TARBALL = cmake-$(CMAKE_VERSION).tar.gz
CMAKE_TARBALL_URL = http://www.cmake.org/files/v2.8/$(CMAKE_TARBALL)
CMAKE_TARBALL_SHA256 = 2b59897864d6220ff20aa8eac64cac8994e004898a1c0f899c8cb4d7b7570b46
CMAKE_CUDA_PATCH = cmake-cuda-2.8.8-0-gc6154b1.patch
CMAKE_CUDA_PATCH_URL = http://sourceforge.net/projects/halmd/files/libs/cmake/$(CMAKE_CUDA_PATCH)
CMAKE_CUDA_PATCH_SHA256 = e2dd8134787543f97ef9399755f4636dc5731ea7c65ef0f7a3bbb1ca48c88b7f
CMAKE_BUILD_DIR = cmake-$(CMAKE_VERSION)
CMAKE_INSTALL_DIR = $(PREFIX)/cmake-$(CMAKE_VERSION)

.fetch-cmake:
	@$(RM) $(CMAKE_TARBALL)
	@$(RM) $(CMAKE_CUDA_PATCH)
	$(WGET) $(CMAKE_TARBALL_URL)
	$(WGET) $(CMAKE_CUDA_PATCH_URL)
	@echo '$(CMAKE_TARBALL_SHA256)  $(CMAKE_TARBALL)' | $(SHA256SUM)
	@echo '$(CMAKE_CUDA_PATCH_SHA256)  $(CMAKE_CUDA_PATCH)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-cmake: .fetch-cmake

.extract-cmake: .fetch-cmake
	$(RM) $(CMAKE_BUILD_DIR)
	$(TAR) -xzf $(CMAKE_TARBALL)
	cd $(CMAKE_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(CMAKE_CUDA_PATCH)
	@$(TOUCH) $@

extract-cmake: .extract-cmake

.configure-cmake: .extract-cmake
	cd $(CMAKE_BUILD_DIR) && ./configure --prefix=$(CMAKE_INSTALL_DIR)
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
	$(RM) $(CMAKE_CUDA_PATCH)

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
LUA_VERSION = 5.2.0
LUA_TARBALL_SHA256 = cabe379465aa8e388988073d59b69e76ba0025429d2c1da80821a252cdf6be0d
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

.fetch-lua:
	@$(RM) $(LUA_TARBALL)
	$(WGET) $(LUA_TARBALL_URL)
	@echo '$(LUA_TARBALL_SHA256)  $(LUA_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-lua: .fetch-lua

.extract-lua: .fetch-lua
	$(RM) $(LUA_BUILD_DIR)
	$(TAR) -xzf $(LUA_TARBALL)
	@$(TOUCH) $@

extract-lua: .extract-lua

.build-lua: .extract-lua
	cd $(LUA_BUILD_DIR) && make linux CFLAGS="$(LUA_CFLAGS)" LIBS="$(LUA_LIBS)"
	@$(TOUCH) $@

build-lua: .build-lua

install-lua: .build-lua
	cd $(LUA_BUILD_DIR) && make install INSTALL_TOP=$(LUA_INSTALL_DIR)

clean-lua:
	@$(RM) .build-lua
	@$(RM) .extract-lua
	$(RM) $(LUA_BUILD_DIR)

distclean-lua: clean-lua
	@$(RM) .fetch-lua
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

LUAJIT_VERSION = 2.0.0-beta9
LUAJIT_TARBALL = LuaJIT-$(LUAJIT_VERSION).tar.gz
LUAJIT_TARBALL_URL = http://luajit.org/download/$(LUAJIT_TARBALL)
LUAJIT_TARBALL_SHA256 = da3793b4364a17c3700d39d13eae799b82ff23da1f61631d735de05333f46240
LUAJIT_PATCH = beta9_hotfix1.patch
LUAJIT_PATCH_URL = http://luajit.org/download/$(LUAJIT_PATCH)
LUAJIT_PATCH_SHA256 = 468234a723c3a2bb7fe8caafc3aac0443473df2790b547a166babaf9b58cc671
LUAJIT_BUILD_DIR = LuaJIT-$(LUAJIT_VERSION)
LUAJIT_INSTALL_DIR = $(PREFIX)/luajit-$(LUAJIT_VERSION)
LUAJIT_CFLAGS = -fPIC

.fetch-luajit:
	@$(RM) $(LUAJIT_TARBALL)
	@$(RM) $(LUAJIT_PATCH)
	$(WGET) $(LUAJIT_TARBALL_URL)
	$(WGET) $(LUAJIT_PATCH_URL)
	@echo '$(LUAJIT_TARBALL_SHA256)  $(LUAJIT_TARBALL)' | $(SHA256SUM)
	@echo '$(LUAJIT_PATCH_SHA256)  $(LUAJIT_PATCH)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-luajit: .fetch-luajit

.extract-luajit: .fetch-luajit
	$(RM) $(LUAJIT_BUILD_DIR)
	$(TAR) -xzf $(LUAJIT_TARBALL)
	cd $(LUAJIT_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(LUAJIT_PATCH)
	@$(TOUCH) $@

extract-luajit: .extract-luajit

.build-luajit: .extract-luajit
	cd $(LUAJIT_BUILD_DIR) && make CFLAGS=$(LUAJIT_CFLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-luajit: .build-luajit

install-luajit: .build-luajit
	cd $(LUAJIT_BUILD_DIR) && make install PREFIX=$(LUAJIT_INSTALL_DIR)
	ln -sf luajit-$(LUAJIT_VERSION) $(LUAJIT_INSTALL_DIR)/bin/lua
	ln -sf libluajit-5.1.a $(LUAJIT_INSTALL_DIR)/lib/liblua.a
	ln -sf luajit-2.0 $(LUAJIT_INSTALL_DIR)/include/lua

clean-luajit:
	@$(RM) .build-luajit
	@$(RM) .extract-luajit
	$(RM) $(LUAJIT_BUILD_DIR)

distclean-luajit: clean-luajit
	@$(RM) .fetch-luajit
	$(RM) $(LUAJIT_TARBALL)
	$(RM) $(LUAJIT_PATCH)

env-luajit:
	@echo
	@echo '# add LuaJIT $(LUAJIT_VERSION) to environment'
	@echo 'export PATH="$(LUAJIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(LUAJIT_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LUAJIT_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Boost C++ libraries with Boost.Log
##

BOOST_VERSION = 1.49.0
BOOST_RELEASE = 1_49_0
BOOST_TARBALL = boost_$(BOOST_RELEASE).tar.bz2
BOOST_TARBALL_URL = http://sourceforge.net/projects/boost/files/boost/$(BOOST_VERSION)/$(BOOST_TARBALL)
BOOST_TARBALL_SHA256 = dd748a7f5507a7e7af74f452e1c52a64e651ed1f7263fce438a06641d2180d3c
BOOST_LOG_VERSION = 1.1
BOOST_LOG_TARBALL = boost-log-$(BOOST_LOG_VERSION).zip
BOOST_LOG_TARBALL_URL = http://sourceforge.net/projects/boost-log/files/boost-log-$(BOOST_LOG_VERSION).zip
BOOST_LOG_TARBALL_SHA256 = 4b00e1d302017298284914c6cc9e7fcae0e097c93e632045d6b0fc4bf6266ba7
BOOST_LOG_DIR = boost-log-$(BOOST_LOG_VERSION)
BOOST_ALWAYS_INLINE_PATCH = boost_1_49_0_fix_integer_log2_always_inline.diff
BOOST_ALWAYS_INLINE_PATCH_URL = http://sourceforge.net/projects/halmd/files/libs/boost/$(BOOST_ALWAYS_INLINE_PATCH)
BOOST_ALWAYS_INLINE_PATCH_SHA256 = 3e0c475d9579804984873858339b53825595c022bfdf811b0a97827cbc18acff
BOOST_LEXICAL_CAST_PATCH = boost_1_49_0_fix_lexical_cast_double_iso_printf.diff
BOOST_LEXICAL_CAST_PATCH_URL = http://sourceforge.net/projects/halmd/files/libs/boost/$(BOOST_LEXICAL_CAST_PATCH)
BOOST_LEXICAL_CAST_PATCH_SHA256 = b738c716a22966d746125078bb3f428ffdad45ae50ee23e434781d904572712c
BOOST_LOG_CXX11_PATCH = boost_log_1.1_fix_compilation_with_c++11.diff
BOOST_LOG_CXX11_PATCH_URL = http://sourceforge.net/projects/halmd/files/libs/boost/$(BOOST_LOG_CXX11_PATCH)
BOOST_LOG_CXX11_PATCH_SHA256 = 350a137852592c5ba079c5fa0a5a5f9c7c2d662fc755d6837c32354c5670ecfe
BOOST_BUILD_DIR = boost_$(BOOST_RELEASE)
BOOST_INSTALL_DIR = $(PREFIX)/boost_$(BOOST_RELEASE)
BOOST_BUILD_FLAGS = cxxflags=-fPIC dll-path=$(BOOST_INSTALL_DIR)/lib

ifndef USE_BZIP2
BOOST_BUILD_FLAGS += -sNO_BZIP2=1
endif
ifndef USE_PYTHON
BOOST_BUILD_FLAGS += --without-python
endif

.fetch-boost:
	@$(RM) $(BOOST_TARBALL)
	@$(RM) $(BOOST_LOG_TARBALL)
	@$(RM) $(BOOST_ALWAYS_INLINE_PATCH)
	@$(RM) $(BOOST_LEXICAL_CAST_PATCH)
	@$(RM) $(BOOST_LOG_CXX11_PATCH)
	$(WGET) $(BOOST_TARBALL_URL)
	$(WGET) $(BOOST_LOG_TARBALL_URL)
	$(WGET) $(BOOST_ALWAYS_INLINE_PATCH_URL)
	$(WGET) $(BOOST_LEXICAL_CAST_PATCH_URL)
	$(WGET) $(BOOST_LOG_CXX11_PATCH_URL)
	@echo '$(BOOST_TARBALL_SHA256)  $(BOOST_TARBALL)' | $(SHA256SUM)
	@echo '$(BOOST_LOG_TARBALL_SHA256)  $(BOOST_LOG_TARBALL)' | $(SHA256SUM)
	@echo '$(BOOST_ALWAYS_INLINE_PATCH_SHA256)  $(BOOST_ALWAYS_INLINE_PATCH)' | $(SHA256SUM)
	@echo '$(BOOST_LEXICAL_CAST_PATCH_SHA256)  $(BOOST_LEXICAL_CAST_PATCH)' | $(SHA256SUM)
	@echo '$(BOOST_LOG_CXX11_PATCH_SHA256)  $(BOOST_LOG_CXX11_PATCH)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-boost: .fetch-boost

.extract-boost: .fetch-boost
	$(RM) $(BOOST_BUILD_DIR) $(BOOST_LOG_DIR)
	$(TAR) -xjf $(BOOST_TARBALL)
	$(UNZIP) $(BOOST_LOG_TARBALL)
	$(CP) $(BOOST_LOG_DIR)/boost/log $(BOOST_BUILD_DIR)/boost/
	$(CP) $(BOOST_LOG_DIR)/libs/log $(BOOST_BUILD_DIR)/libs/
	cd $(BOOST_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(BOOST_ALWAYS_INLINE_PATCH)
	cd $(BOOST_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(BOOST_LEXICAL_CAST_PATCH)
	cd $(BOOST_BUILD_DIR) && $(PATCH) -p3 < $(CURDIR)/$(BOOST_LOG_CXX11_PATCH)
	@$(TOUCH) $@

extract-boost: .extract-boost

.configure-boost: .extract-boost
	cd $(BOOST_BUILD_DIR) && ./bootstrap.sh
	@$(TOUCH) $@

configure-boost: .configure-boost

.build-boost: .configure-boost
	cd $(BOOST_BUILD_DIR) && ./bjam $(BOOST_BUILD_FLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-boost: .build-boost

install-boost: .build-boost
	cd $(BOOST_BUILD_DIR) && ./bjam $(BOOST_BUILD_FLAGS) install --prefix=$(BOOST_INSTALL_DIR)

clean-boost:
	@$(RM) .build-boost
	@$(RM) .configure-boost
	@$(RM) .extract-boost
	$(RM) $(BOOST_BUILD_DIR) $(BOOST_LOG_DIR)

distclean-boost: clean-boost
	@$(RM) .fetch-boost
	$(RM) $(BOOST_TARBALL)
	$(RM) $(BOOST_LOG_TARBALL)
	$(RM) $(BOOST_BOOST_ALWAYS_INLINE_PATCH)
	$(RM) $(BOOST_BOOST_LEXICAL_CAST_PATCH)
	$(RM) $(BOOST_BOOST_LOG_CXX11_PATCH)

env-boost:
	@echo
	@echo '# add Boost $(BOOST_VERSION) to environment'
	@echo 'export PATH="$(BOOST_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(BOOST_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export PYTHONPATH="$(BOOST_INSTALL_DIR)/lib$${PYTHONPATH+:$$PYTHONPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(BOOST_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Luabind library
##

LUABIND_VERSION = 0.9.1-29-g5778ff8
LUABIND_TARBALL = luabind-$(LUABIND_VERSION).tar.bz2
LUABIND_TARBALL_URL = http://sourceforge.net/projects/halmd/files/libs/luabind/$(LUABIND_TARBALL)
LUABIND_TARBALL_SHA256 = 9d85fa51e9693a8050d0cd9e0eaceb045db1a86f38209102bb0c120d88f318e6
LUABIND_BUILD_DIR = luabind-$(LUABIND_VERSION)
LUABIND_BUILD_FLAGS = cxxflags=-fPIC link=static variant=release variant=debug
LUABIND_INSTALL_DIR = $(PREFIX)/luabind-$(LUABIND_VERSION)

.fetch-luabind:
	@$(RM) $(LUABIND_TARBALL)
	$(WGET) $(LUABIND_TARBALL_URL)
	@echo '$(LUABIND_TARBALL_SHA256)  $(LUABIND_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-luabind: .fetch-luabind

ifdef USE_LUAJIT
LUABIND_BUILD_LUA_TARGET = .build-luajit
else
LUABIND_BUILD_LUA_TARGET = .build-lua
endif

.extract-luabind: .fetch-luabind $(LUABIND_BUILD_LUA_TARGET)
	$(RM) $(LUABIND_BUILD_DIR)
	$(TAR) -xjf $(LUABIND_TARBALL)
	mkdir $(LUABIND_BUILD_DIR)/lua
ifdef USE_LUAJIT
	mkdir $(LUABIND_BUILD_DIR)/lua/lib
	ln -s $(CURDIR)/$(LUAJIT_BUILD_DIR)/src $(LUABIND_BUILD_DIR)/lua/include
	ln -s $(CURDIR)/$(LUAJIT_BUILD_DIR)/src/libluajit.a $(LUABIND_BUILD_DIR)/lua/lib/liblua.a
else
	ln -s $(CURDIR)/$(LUA_BUILD_DIR)/src $(LUABIND_BUILD_DIR)/lua/include
	ln -s $(CURDIR)/$(LUA_BUILD_DIR)/src $(LUABIND_BUILD_DIR)/lua/lib
endif
	@$(TOUCH) $@

extract-luabind: .extract-luabind

.build-luabind: .extract-luabind .configure-boost
	cd $(LUABIND_BUILD_DIR) && BOOST_ROOT=$(CURDIR)/$(BOOST_BUILD_DIR) LUA_PATH=$(CURDIR)/$(LUABIND_BUILD_DIR)/lua $(CURDIR)/$(BOOST_BUILD_DIR)/bjam $(LUABIND_BUILD_FLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-luabind: .build-luabind

install-luabind: .build-luabind
	cd $(LUABIND_BUILD_DIR) && BOOST_ROOT=$(CURDIR)/$(BOOST_BUILD_DIR) LUA_PATH=$(CURDIR)/$(LUABIND_BUILD_DIR)/lua $(CURDIR)/$(BOOST_BUILD_DIR)/bjam $(LUABIND_BUILD_FLAGS) install --prefix=$(LUABIND_INSTALL_DIR)

clean-luabind:
	@$(RM) .build-luabind
	@$(RM) .extract-luabind
	$(RM) $(LUABIND_BUILD_DIR)

distclean-luabind: clean-luabind
	@$(RM) .fetch-luabind
	$(RM) $(LUABIND_TARBALL)

env-luabind:
	@echo
	@echo '# add Luabind $(LUABIND_VERSION) to environment'
	@echo 'export LD_LIBRARY_PATH="$(LUABIND_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(LUABIND_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## HDF5 C++ library
##

HDF5_VERSION = 1.8.8
HDF5_TARBALL = hdf5-$(HDF5_VERSION).tar.bz2
HDF5_TARBALL_URL = http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_VERSION)/src/$(HDF5_TARBALL)
HDF5_TARBALL_SHA256 = b0ebb0b5478c6c0427631d4ad08f96e39f1b09fde615aa98d2a1b8fb7f6dced3
HDF5_BUILD_DIR = hdf5-$(HDF5_VERSION)
HDF5_INSTALL_DIR = $(PREFIX)/hdf5-$(HDF5_VERSION)
HDF5_CONFIGURE_FLAGS = --enable-cxx
HDF5_CFLAGS = -fPIC
HDF5_CXXFLAGS = -fPIC

.fetch-hdf5:
	@$(RM) $(HDF5_TARBALL)
	$(WGET) $(HDF5_TARBALL_URL)
	@echo '$(HDF5_TARBALL_SHA256)  $(HDF5_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-hdf5: .fetch-hdf5

.extract-hdf5: .fetch-hdf5
	$(RM) $(HDF5_BUILD_DIR)
	$(TAR) -xjf $(HDF5_TARBALL)
	@$(TOUCH) $@

extract-hdf5: .extract-hdf5

.configure-hdf5: .extract-hdf5
	cd $(HDF5_BUILD_DIR) && CFLAGS="$(HDF5_CFLAGS)" CXXFLAGS="$(HDF5_CXXFLAGS)" ./configure $(HDF5_CONFIGURE_FLAGS) --prefix=$(HDF5_INSTALL_DIR)
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
	@echo 'export PATH="$(HDF5_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(HDF5_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(HDF5_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## Git version control
##

GIT_VERSION = 1.7.9.1
GIT_TARBALL = git-$(GIT_VERSION).tar.gz
GIT_TARBALL_URL = http://git-core.googlecode.com/files/$(GIT_TARBALL)
GIT_TARBALL_SHA256 = b689a0ddbc99f9a69aef7c81c569289a28ba0787cd27e5e188112e1c3f0e8152
GIT_MANPAGES_TARBALL = git-manpages-$(GIT_VERSION).tar.gz
GIT_MANPAGES_TARBALL_URL = http://git-core.googlecode.com/files/$(GIT_MANPAGES_TARBALL)
GIT_MANPAGES_TARBALL_SHA256 = f31d91061e96b5f882ceed2160d44937b2679931e7e217a66bfe3a23df46adae
GIT_BUILD_DIR = git-$(GIT_VERSION)
GIT_CONFIGURE_FLAGS = --without-python
GIT_INSTALL_DIR = $(PREFIX)/git-$(GIT_VERSION)

.fetch-git:
	@$(RM) $(GIT_TARBALL)
	@$(RM) $(GIT_MANPAGES_TARBALL)
	$(WGET) $(GIT_TARBALL_URL)
	$(WGET) $(GIT_MANPAGES_TARBALL_URL)
	@echo '$(GIT_TARBALL_SHA256)  $(GIT_TARBALL)' | $(SHA256SUM)
	@echo '$(GIT_MANPAGES_TARBALL_SHA256)  $(GIT_MANPAGES_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-git: .fetch-git

.extract-git: .fetch-git
	$(RM) $(GIT_BUILD_DIR)
	$(TAR) -xzf $(GIT_TARBALL)
	@$(TOUCH) $@

extract-git: .extract-git

.configure-git: .extract-git
	cd $(GIT_BUILD_DIR) && ./configure $(GIT_CONFIGURE_FLAGS) --prefix=$(GIT_INSTALL_DIR)
	@$(TOUCH) $@

configure-git: .configure-git

.build-git: .configure-git
	cd $(GIT_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-git: .build-git

install-git: .build-git
	cd $(GIT_BUILD_DIR) && make install
	install -d $(GIT_INSTALL_DIR)/share/man
	cd $(GIT_INSTALL_DIR)/share/man && $(TAR) -xzf $(CURDIR)/$(GIT_MANPAGES_TARBALL)

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
	@echo 'export PATH="$(GIT_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GIT_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## htop
##

HTOP_VERSION = 1.0.1
HTOP_TARBALL = htop-$(HTOP_VERSION).tar.gz
HTOP_TARBALL_URL = http://sourceforge.net/projects/htop/files/htop/$(HTOP_VERSION)/$(HTOP_TARBALL)
HTOP_TARBALL_SHA256 = 07db2cbe02835f9e186b9610ecc3beca330a5c9beadb3b6069dd0a10561506f2
HTOP_BUILD_DIR = htop-$(HTOP_VERSION)
HTOP_INSTALL_DIR = $(PREFIX)/htop-$(HTOP_VERSION)

.fetch-htop:
	@$(RM) $(HTOP_TARBALL)
	$(WGET) $(HTOP_TARBALL_URL)
	@echo '$(HTOP_TARBALL_SHA256)  $(HTOP_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-htop: .fetch-htop

.extract-htop: .fetch-htop
	$(RM) $(HTOP_BUILD_DIR)
	$(TAR) -xzf $(HTOP_TARBALL)
	@$(TOUCH) $@

extract-htop: .extract-htop

.configure-htop: .extract-htop
	cd $(HTOP_BUILD_DIR) && ./configure --prefix=$(HTOP_INSTALL_DIR)
	@$(TOUCH) $@

configure-htop: .configure-htop

.build-htop: .configure-htop
	cd $(HTOP_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-htop: .build-htop

install-htop: .build-htop
	cd $(HTOP_BUILD_DIR) && make install

clean-htop:
	@$(RM) .build-htop
	@$(RM) .configure-htop
	@$(RM) .extract-htop
	$(RM) $(HTOP_BUILD_DIR)

distclean-htop: clean-htop
	@$(RM) .fetch-htop
	$(RM) $(HTOP_TARBALL)

env-htop:
	@echo
	@echo '# add htop $(HTOP_VERSION) to environment'
	@echo 'export PATH="$(HTOP_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(HTOP_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## python-sphinx
##

PYTHON_SPHINX_VERSION = 1.1.2
PYTHON_SPHINX_TARBALL = Sphinx-$(PYTHON_SPHINX_VERSION).tar.gz
PYTHON_SPHINX_TARBALL_URL = http://pypi.python.org/packages/source/S/Sphinx/$(PYTHON_SPHINX_TARBALL)
PYTHON_SPHINX_TARBALL_SHA256 = cf66ee61eef61b7c478907282bddcdd5e04eebd69a00a2bb93881427938fe688
PYTHON_SPHINX_BUILD_DIR = Sphinx-$(PYTHON_SPHINX_VERSION)
PYTHON_SPHINX_INSTALL_DIR = $(PREFIX)/python-sphinx-$(PYTHON_SPHINX_VERSION)
PYTHON_SPHINX_PYTHONPATH = $(PYTHON_SPHINX_INSTALL_DIR)/lib/python

.fetch-python-sphinx:
	@$(RM) $(PYTHON_SPHINX_TARBALL)
	$(WGET) $(PYTHON_SPHINX_TARBALL_URL)
	@echo '$(PYTHON_SPHINX_TARBALL_SHA256)  $(PYTHON_SPHINX_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-python-sphinx: .fetch-python-sphinx

.extract-python-sphinx: .fetch-python-sphinx
	$(RM) $(PYTHON_SPHINX_BUILD_DIR)
	$(TAR) -xzf $(PYTHON_SPHINX_TARBALL)
	@$(TOUCH) $@

extract-python-sphinx: .extract-python-sphinx

.build-python-sphinx: .extract-python-sphinx
	cd $(PYTHON_SPHINX_BUILD_DIR) && python setup.py build
	@$(TOUCH) $@

build-python-sphinx: .build-python-sphinx

install-python-sphinx: .build-python-sphinx
	install -d $(PYTHON_SPHINX_PYTHONPATH)
	cd $(PYTHON_SPHINX_BUILD_DIR) && PYTHONPATH="$(PYTHON_SPHINX_PYTHONPATH)$${PYTHONPATH+:$$PYTHONPATH}" python setup.py install --home=$(PYTHON_SPHINX_INSTALL_DIR)

clean-python-sphinx:
	@$(RM) .build-python-sphinx
	@$(RM) .extract-python-sphinx
	$(RM) $(PYTHON_SPHINX_BUILD_DIR)

distclean-python-sphinx: clean-python-sphinx
	@$(RM) .fetch-python-sphinx
	$(RM) $(PYTHON_SPHINX_TARBALL)

env-python-sphinx:
	@echo
	@echo '# add python-sphinx $(PYTHON_SPHINX_VERSION) to environment'
	@echo 'export PATH="$(PYTHON_SPHINX_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export PYTHONPATH="$(PYTHON_SPHINX_PYTHONPATH)$${PYTHONPATH+:$$PYTHONPATH}"'

##
## Doxygen
##

DOXYGEN_VERSION = 1.7.6.1
DOXYGEN_TARBALL = doxygen-$(DOXYGEN_VERSION).src.tar.gz
DOXYGEN_TARBALL_URL = http://ftp.stack.nl/pub/users/dimitri/$(DOXYGEN_TARBALL)
DOXYGEN_TARBALL_SHA256 = 0e60e794fb172d3fa4a9a9535f0b8e0eeb04e8366153f6b417569af0bcd61fcd
DOXYGEN_BUILD_DIR = doxygen-$(DOXYGEN_VERSION)
DOXYGEN_INSTALL_DIR = $(PREFIX)/doxygen-$(DOXYGEN_VERSION)

.fetch-doxygen:
	@$(RM) $(DOXYGEN_TARBALL)
	$(WGET) $(DOXYGEN_TARBALL_URL)
	@echo '$(DOXYGEN_TARBALL_SHA256)  $(DOXYGEN_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-doxygen: .fetch-doxygen

.extract-doxygen: .fetch-doxygen
	$(RM) $(DOXYGEN_BUILD_DIR)
	$(TAR) -xzf $(DOXYGEN_TARBALL)
	@$(TOUCH) $@

extract-doxygen: .extract-doxygen

.configure-doxygen: .extract-doxygen
	cd $(DOXYGEN_BUILD_DIR) && ./configure --prefix $(DOXYGEN_INSTALL_DIR)
	@$(TOUCH) $@

configure-doxygen: .configure-doxygen

.build-doxygen: .configure-doxygen
	cd $(DOXYGEN_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-doxygen: .build-doxygen

install-doxygen: .build-doxygen
	cd $(DOXYGEN_BUILD_DIR) && make install

clean-doxygen:
	@$(RM) .build-doxygen
	@$(RM) .configure-doxygen
	@$(RM) .extract-doxygen
	$(RM) $(DOXYGEN_BUILD_DIR)

distclean-doxygen: clean-doxygen
	@$(RM) .fetch-doxygen
	$(RM) $(DOXYGEN_TARBALL)

env-doxygen:
	@echo
	@echo '# add Doxygen $(DOXYGEN_VERSION) to environment'
	@echo 'export PATH="$(DOXYGEN_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(DOXYGEN_INSTALL_DIR)/man$${MANPATH+:$$MANPATH}"'

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

.fetch-graphviz:
	@$(RM) $(GRAPHVIZ_TARBALL)
	$(WGET) $(GRAPHVIZ_TARBALL_URL)
	@echo '$(GRAPHVIZ_TARBALL_SHA256)  $(GRAPHVIZ_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-graphviz: .fetch-graphviz

.extract-graphviz: .fetch-graphviz
	$(RM) $(GRAPHVIZ_BUILD_DIR)
	$(TAR) -xzf $(GRAPHVIZ_TARBALL)
	@$(TOUCH) $@

extract-graphviz: .extract-graphviz

.configure-graphviz: .extract-graphviz
	cd $(GRAPHVIZ_BUILD_DIR) && ./configure $(GRAPHVIZ_CONFIGURE_FLAGS) --prefix=$(GRAPHVIZ_INSTALL_DIR)
	@$(TOUCH) $@

configure-graphviz: .configure-graphviz

.build-graphviz: .configure-graphviz
	cd $(GRAPHVIZ_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-graphviz: .build-graphviz

install-graphviz: .build-graphviz
	cd $(GRAPHVIZ_BUILD_DIR) && make install

clean-graphviz:
	@$(RM) .build-graphviz
	@$(RM) .configure-graphviz
	@$(RM) .extract-graphviz
	$(RM) $(GRAPHVIZ_BUILD_DIR)

distclean-graphviz: clean-graphviz
	@$(RM) .fetch-graphviz
	$(RM) $(GRAPHVIZ_TARBALL)

env-graphviz:
	@echo
	@echo '# add Graphviz $(GRAPHVIZ_VERSION) to environment'
	@echo 'export PATH="$(GRAPHVIZ_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GRAPHVIZ_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## Clang C++ compiler
##

CLANG_VERSION = 3.0
LLVM_TARBALL = llvm-$(CLANG_VERSION).tar.gz
LLVM_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(LLVM_TARBALL)
LLVM_TARBALL_SHA256 = 519eb11d3499ce99c6ffdb8718651fc91425ed7690eac91c8d6853474f7c0477
CLANG_TARBALL = clang-$(CLANG_VERSION).tar.gz
CLANG_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(CLANG_TARBALL)
CLANG_TARBALL_SHA256 = b64e72da356d7c3428cfd7ac620d49ec042c84eaee13c26024879f555f4e126d
CLANG_BUILD_DIR = llvm-$(CLANG_VERSION).src
CLANG_CONFIGURE_FLAGS = --enable-optimized --enable-bindings=none
CLANG_INSTALL_DIR = $(PREFIX)/clang-$(CLANG_VERSION)

.fetch-clang:
	@$(RM) $(LLVM_TARBALL)
	@$(RM) $(CLANG_TARBALL)
	$(WGET) $(LLVM_TARBALL_URL)
	$(WGET) $(CLANG_TARBALL_URL)
	@echo '$(LLVM_TARBALL_SHA256)  $(LLVM_TARBALL)' | $(SHA256SUM)
	@echo '$(CLANG_TARBALL_SHA256)  $(CLANG_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-clang: .fetch-clang

.extract-clang: .fetch-clang
	$(RM) $(CLANG_BUILD_DIR)
	$(TAR) -xzf $(LLVM_TARBALL)
	cd $(CLANG_BUILD_DIR)/tools && $(TAR) -xzf $(CURDIR)/$(CLANG_TARBALL) && mv clang-$(CLANG_VERSION).src clang
	@$(TOUCH) $@

extract-clang: .extract-clang

.configure-clang: .extract-clang
	cd $(CLANG_BUILD_DIR) && ./configure $(CLANG_CONFIGURE_FLAGS) --prefix=$(CLANG_INSTALL_DIR)
	@$(TOUCH) $@

configure-clang: .configure-clang

.build-clang: .configure-clang
	cd $(CLANG_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-clang: .build-clang

install-clang: .build-clang
	cd $(CLANG_BUILD_DIR) && make install

clean-clang:
	@$(RM) .build-clang
	@$(RM) .configure-clang
	@$(RM) .extract-clang
	$(RM) $(CLANG_BUILD_DIR)

distclean-clang: clean-clang
	@$(RM) .fetch-clang
	$(RM) $(LLVM_TARBALL)
	$(RM) $(CLANG_TARBALL)

env-clang:
	@echo
	@echo '# add Clang $(CLANG_VERSION) to environment'
	@echo 'export PATH="$(CLANG_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(CLANG_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## GNU Parallel
##

GNU_PARALLEL_VERSION = 20111222
GNU_PARALLEL_TARBALL = parallel-$(GNU_PARALLEL_VERSION).tar.bz2
GNU_PARALLEL_TARBALL_URL = http://ftp.gnu.org/gnu/parallel/$(GNU_PARALLEL_TARBALL)
GNU_PARALLEL_TARBALL_SHA256 = 2bb8eb1f358963eb50d3f9f285b887c378e0a660061e7e3e9b8205bf2c443766
GNU_PARALLEL_BUILD_DIR = parallel-$(GNU_PARALLEL_VERSION)
GNU_PARALLEL_INSTALL_DIR = $(PREFIX)/parallel-$(GNU_PARALLEL_VERSION)

.fetch-gnu-parallel:
	@$(RM) $(GNU_PARALLEL_TARBALL)
	$(WGET) $(GNU_PARALLEL_TARBALL_URL)
	@echo '$(GNU_PARALLEL_TARBALL_SHA256)  $(GNU_PARALLEL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gnu-parallel: .fetch-gnu-parallel

.extract-gnu-parallel: .fetch-gnu-parallel
	$(RM) $(GNU_PARALLEL_BUILD_DIR)
	$(TAR) -xjf $(GNU_PARALLEL_TARBALL)
	@$(TOUCH) $@

extract-gnu-parallel: .extract-gnu-parallel

.configure-gnu-parallel: .extract-gnu-parallel
	cd $(GNU_PARALLEL_BUILD_DIR) && ./configure --prefix=$(GNU_PARALLEL_INSTALL_DIR)
	@$(TOUCH) $@

configure-gnu-parallel: .configure-gnu-parallel

.build-gnu-parallel: .configure-gnu-parallel
	cd $(GNU_PARALLEL_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gnu-parallel: .build-gnu-parallel

install-gnu-parallel: .build-gnu-parallel
	cd $(GNU_PARALLEL_BUILD_DIR) && make install

clean-gnu-parallel:
	@$(RM) .build-gnu-parallel
	@$(RM) .configure-gnu-parallel
	@$(RM) .extract-gnu-parallel
	$(RM) $(GNU_PARALLEL_BUILD_DIR)

distclean-gnu-parallel: clean-gnu-parallel
	@$(RM) .fetch-gnu-parallel
	$(RM) $(GNU_PARALLEL_TARBALL)

env-gnu-parallel:
	@echo
	@echo '# add GNU Parallel $(GNU_PARALLEL_VERSION) to environment'
	@echo 'export PATH="$(GNU_PARALLEL_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GNU_PARALLEL_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## GCC (GNU Compiler Collection)
##

GCC_VERSION = 4.7.0
GCC_TARBALL = gcc-$(GCC_VERSION).tar.bz2
GCC_TARBALL_URL = http://mirror.csclub.uwaterloo.ca/gnu/gcc/gcc-$(GCC_VERSION)/$(GCC_TARBALL)
GCC_TARBALL_SHA256 = a680083e016f656dab7acd45b9729912e70e71bbffcbf0e3e8aa1cccf19dc9a5
GCC_BUILD_DIR = gcc-$(GCC_VERSION)
GCC_INSTALL_DIR = $(PREFIX)/gcc-$(GCC_VERSION)

.fetch-gcc:
	@$(RM) $(GCC_TARBALL)
	$(WGET) $(GCC_TARBALL_URL)
	@echo '$(GCC_TARBALL_SHA256)  $(GCC_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gcc: .fetch-gcc .fetch-gmp .fetch-mpfr .fetch-mpc .fetch-ppl .fetch-cloog-ppl

.extract-gcc: .fetch-gcc
	$(RM) $(GCC_BUILD_DIR)
	$(TAR) -xjf $(GCC_TARBALL)
	@$(TOUCH) $@

extract-gcc: .extract-gcc

.configure-gcc: .extract-gcc .install-gmp .install-mpfr .install-mpc .install-ppl .install-cloog-ppl
	cd $(GCC_BUILD_DIR) && LDFLAGS=-L$(GMP_INSTALL_DIR)/lib ac_cv_lib_pwl_PWL_handle_timeout=no ./configure --prefix=$(GCC_INSTALL_DIR) --disable-multilib --enable-languages=c,c++,fortran,lto --with-gmp=$(GMP_INSTALL_DIR) --with-mpfr=$(MPFR_INSTALL_DIR) --with-mpc=$(MPC_INSTALL_DIR) --with-ppl=$(PPL_INSTALL_DIR) --with-cloog=$(CLOOG_PPL_INSTALL_DIR) --enable-build-with-cxx
	@$(TOUCH) $@

configure-gcc: .configure-gcc

.build-gcc: .configure-gcc
	cd $(GCC_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gcc: .build-gcc

install-gcc: .build-gcc
	cd $(GCC_BUILD_DIR) && make install

clean-gcc:
	@$(RM) .build-gcc
	@$(RM) .configure-gcc
	@$(RM) .extract-gcc
	@$(RM) .install-gmp
	@$(RM) .install-mpfr
	@$(RM) .install-mpc
	@$(RM) .install-ppl
	@$(RM) .install-cloog-ppl
	$(RM) $(GCC_BUILD_DIR)

distclean-gcc: clean-gcc
	@$(RM) .fetch-gcc
	$(RM) $(GCC_TARBALL)

env-gcc:
	@echo
	@echo '# add gcc $(GCC_VERSION) to environment'
	@echo 'export PATH="$(GCC_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GCC_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'
	@echo 'export LD_LIBRARY_PATH="$(GCC_INSTALL_DIR)/lib64$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'

##
## GMP (GNU Multiple Precision Arithmetic Library)
##

GMP_VERSION = 5.0.4
GMP_TARBALL = gmp-$(GMP_VERSION).tar.bz2
GMP_TARBALL_URL = ftp://ftp.gmplib.org/pub/gmp-$(GMP_VERSION)/$(GMP_TARBALL)
GMP_TARBALL_SHA256 = 35d4aade3e4bdf0915c944599b10d23f108ffedf6c3188aeec52221c5cf9a06f
GMP_BUILD_DIR = gmp-$(GMP_VERSION)
GMP_INSTALL_DIR = $(CURDIR)/$(GCC_BUILD_DIR)/gmp

.fetch-gmp:
	@$(RM) $(GMP_TARBALL)
	$(WGET) $(GMP_TARBALL_URL)
	@echo '$(GMP_TARBALL_SHA256)  $(GMP_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-gmp: .fetch-gmp

.extract-gmp: .fetch-gmp
	$(RM) $(GMP_BUILD_DIR)
	$(TAR) -xjf $(GMP_TARBALL)
	@$(TOUCH) $@

extract-gmp: .extract-gmp

.configure-gmp: .extract-gmp
	cd $(GMP_BUILD_DIR) && ./configure --prefix=$(GMP_INSTALL_DIR) --enable-cxx --disable-shared
	@$(TOUCH) $@

configure-gmp: .configure-gmp

.build-gmp: .configure-gmp
	cd $(GMP_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gmp: .build-gmp

.install-gmp: .build-gmp
	cd $(GMP_BUILD_DIR) && make install
	@$(TOUCH) $@

clean-gmp:
	@$(RM) .build-gmp
	@$(RM) .configure-gmp
	@$(RM) .extract-gmp
	$(RM) $(GMP_BUILD_DIR)

distclean-gmp: clean-gmp
	@$(RM) .fetch-gmp
	$(RM) $(GMP_TARBALL)

##
## MPFR (Multiple-precision floating-point computations with correct rounding)
##

MPFR_VERSION = 3.1.0
MPFR_TARBALL = mpfr-$(MPFR_VERSION).tar.bz2
MPFR_TARBALL_URL = http://www.mpfr.org/mpfr-$(MPFR_VERSION)/$(MPFR_TARBALL)
MPFR_TARBALL_SHA256 = 74a7bbbad168dd1cc414f1c9210b8fc16ccfc8e422d34b3371a8978e31eab680
MPFR_BUILD_DIR = mpfr-$(MPFR_VERSION)
MPFR_INSTALL_DIR = $(CURDIR)/$(GCC_BUILD_DIR)/mpfr

.fetch-mpfr:
	@$(RM) $(MPFR_TARBALL)
	$(WGET) $(MPFR_TARBALL_URL)
	@echo '$(MPFR_TARBALL_SHA256)  $(MPFR_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-mpfr: .fetch-mpfr

.extract-mpfr: .fetch-mpfr
	$(RM) $(MPFR_BUILD_DIR)
	$(TAR) -xjf $(MPFR_TARBALL)
	@$(TOUCH) $@

extract-mpfr: .extract-mpfr

.configure-mpfr: .extract-mpfr .install-gmp
	cd $(MPFR_BUILD_DIR) && ./configure --prefix=$(MPFR_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --disable-shared
	@$(TOUCH) $@

configure-mpfr: .configure-mpfr

.build-mpfr: .configure-mpfr
	cd $(MPFR_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-mpfr: .build-mpfr

.install-mpfr: .build-mpfr
	cd $(MPFR_BUILD_DIR) && make install
	@$(TOUCH) $@

clean-mpfr:
	@$(RM) .build-mpfr
	@$(RM) .configure-mpfr
	@$(RM) .extract-mpfr
	$(RM) $(MPFR_BUILD_DIR)

distclean-mpfr: clean-mpfr
	@$(RM) .fetch-mpfr
	$(RM) $(MPFR_TARBALL)

##
## MPC (arithmetic of complex numbers with arbitrarily high precision and correct rounding)
##

MPC_VERSION = 0.9
MPC_TARBALL = mpc-$(MPC_VERSION).tar.gz
MPC_TARBALL_URL = http://www.multiprecision.org/mpc/download/$(MPC_TARBALL)
MPC_TARBALL_SHA256 = fd3efe422f0d454592059e80f2c00d1a2e381bf2beda424c5094abd4deb049ac
MPC_BUILD_DIR = mpc-$(MPC_VERSION)
MPC_INSTALL_DIR = $(CURDIR)/$(GCC_BUILD_DIR)/mpc

.fetch-mpc:
	@$(RM) $(MPC_TARBALL)
	$(WGET) $(MPC_TARBALL_URL)
	@echo '$(MPC_TARBALL_SHA256)  $(MPC_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-mpc: .fetch-mpc

.extract-mpc: .fetch-mpc
	$(RM) $(MPC_BUILD_DIR)
	$(TAR) -xzf $(MPC_TARBALL)
	@$(TOUCH) $@

extract-mpc: .extract-mpc

.configure-mpc: .extract-mpc .install-gmp .install-mpfr
	cd $(MPC_BUILD_DIR) && ./configure --prefix=$(MPC_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --with-mpfr=$(MPFR_INSTALL_DIR) --disable-shared
	@$(TOUCH) $@

configure-mpc: .configure-mpc

.build-mpc: .configure-mpc
	cd $(MPC_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-mpc: .build-mpc

.install-mpc: .build-mpc
	cd $(MPC_BUILD_DIR) && make install
	@$(TOUCH) $@

clean-mpc:
	@$(RM) .build-mpc
	@$(RM) .configure-mpc
	@$(RM) .extract-mpc
	$(RM) $(MPC_BUILD_DIR)

distclean-mpc: clean-mpc
	@$(RM) .fetch-mpc
	$(RM) $(MPC_TARBALL)

##
## PPL (Parma Polyhedra Library)
##

PPL_VERSION = 0.11.2
PPL_TARBALL = ppl-$(PPL_VERSION).tar.bz2
PPL_TARBALL_URL = ftp://ftp.cs.unipr.it/pub/ppl/releases/$(PPL_VERSION)/$(PPL_TARBALL)
PPL_TARBALL_SHA256 = e3fbd1c19ef44c6f020951807cdb6fc6a8153cd3a5c53b0ab9cf4c4f6e8cbbeb
PPL_BUILD_DIR = ppl-$(PPL_VERSION)
PPL_INSTALL_DIR = $(CURDIR)/$(GCC_BUILD_DIR)/ppl

.fetch-ppl:
	@$(RM) $(PPL_TARBALL)
	$(WGET) $(PPL_TARBALL_URL)
	@echo '$(PPL_TARBALL_SHA256)  $(PPL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-ppl: .fetch-ppl

.extract-ppl: .fetch-ppl
	$(RM) $(PPL_BUILD_DIR)
	$(TAR) -xjf $(PPL_TARBALL)
	@$(TOUCH) $@

extract-ppl: .extract-ppl

.configure-ppl: .extract-ppl .install-gmp
	cd $(PPL_BUILD_DIR) && CPPFLAGS=-I$(GMP_INSTALL_DIR)/include LDFLAGS=-L$(GMP_INSTALL_DIR)/lib ./configure --prefix=$(PPL_INSTALL_DIR) --enable-interfaces=c,cxx --disable-shared --disable-watchdog
	@$(TOUCH) $@

configure-ppl: .configure-ppl

.build-ppl: .configure-ppl
	cd $(PPL_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-ppl: .build-ppl

.install-ppl: .build-ppl
	cd $(PPL_BUILD_DIR) && make install
	@$(TOUCH) $@

clean-ppl:
	@$(RM) .build-ppl
	@$(RM) .configure-ppl
	@$(RM) .extract-ppl
	$(RM) $(PPL_BUILD_DIR)

distclean-ppl: clean-ppl
	@$(RM) .fetch-ppl
	$(RM) $(PPL_TARBALL)

##
## CLooG-PPL
##

CLOOG_PPL_VERSION = 0.15.11
CLOOG_PPL_TARBALL = cloog-ppl-$(CLOOG_PPL_VERSION).tar.gz
CLOOG_PPL_TARBALL_URL = ftp://ftp.fu-berlin.de/unix/languages/gcc/infrastructure/$(CLOOG_PPL_TARBALL)
CLOOG_PPL_TARBALL_SHA256 = 7cd634d0b2b401b04096b545915ac67f883556e9a524e8e803a6bf6217a84d5f
CLOOG_PPL_BUILD_DIR = cloog-ppl-$(CLOOG_PPL_VERSION)
CLOOG_PPL_INSTALL_DIR = $(CURDIR)/$(GCC_BUILD_DIR)/cloog-ppl

.fetch-cloog-ppl:
	@$(RM) $(CLOOG_PPL_TARBALL)
	$(WGET) $(CLOOG_PPL_TARBALL_URL)
	@echo '$(CLOOG_PPL_TARBALL_SHA256)  $(CLOOG_PPL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-cloog-ppl: .fetch-cloog-ppl

.extract-cloog-ppl: .fetch-cloog-ppl
	$(RM) $(CLOOG_PPL_BUILD_DIR)
	$(TAR) -xzf $(CLOOG_PPL_TARBALL)
	@$(TOUCH) $@

extract-cloog-ppl: .extract-cloog-ppl

.configure-cloog-ppl: .extract-cloog-ppl .install-ppl .install-gmp
	cd $(CLOOG_PPL_BUILD_DIR) && ./configure --prefix=$(CLOOG_PPL_INSTALL_DIR) --with-ppl=$(PPL_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --disable-shared --with-host-libstdcxx=-lstdc++
	@$(TOUCH) $@

configure-cloog-ppl: .configure-cloog-ppl

.build-cloog-ppl: .configure-cloog-ppl
	cd $(CLOOG_PPL_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-cloog-ppl: .build-cloog-ppl

.install-cloog-ppl: .build-cloog-ppl
	cd $(CLOOG_PPL_BUILD_DIR) && make install
	@$(TOUCH) $@

clean-cloog-ppl:
	@$(RM) .build-cloog-ppl
	@$(RM) .configure-cloog-ppl
	@$(RM) .extract-cloog-ppl
	$(RM) $(CLOOG_PPL_BUILD_DIR)

distclean-cloog-ppl: clean-cloog-ppl
	@$(RM) .fetch-cloog-ppl
	$(RM) $(CLOOG_PPL_TARBALL)
