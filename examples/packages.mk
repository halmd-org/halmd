#
# Copyright © 2011-2012 Peter Colberg
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
GIT = git
TAR = tar
RM = rm -rf
CP = cp -r
TOUCH = touch
PATCH = patch
SHA256SUM = sha256sum --check

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

CMAKE_VERSION = 2.8.8
CMAKE_TARBALL = cmake-$(CMAKE_VERSION).tar.gz
CMAKE_TARBALL_URL = http://www.cmake.org/files/v2.8/$(CMAKE_TARBALL)
CMAKE_TARBALL_SHA256 = 2b59897864d6220ff20aa8eac64cac8994e004898a1c0f899c8cb4d7b7570b46
CMAKE_CUDA_VERSION = ${CMAKE_VERSION}-2-ga06ff77
CMAKE_CUDA_PATCH = cmake-cuda-${CMAKE_CUDA_VERSION}.patch
CMAKE_CUDA_PATCH_URL = http://sourceforge.net/projects/halmd/files/libs/cmake/$(CMAKE_CUDA_PATCH)
CMAKE_CUDA_PATCH_SHA256 = 04b22de5479f8a869b2e4dcee08de69b2083d645ea5c1b35cef366c2d174c941
CMAKE_BUILD_DIR = cmake-$(CMAKE_VERSION)
CMAKE_INSTALL_DIR = $(PREFIX)/cmake-$(CMAKE_CUDA_VERSION)

.fetch-cmake-$(CMAKE_CUDA_VERSION):
	@$(RM) $(CMAKE_TARBALL)
	@$(RM) $(CMAKE_CUDA_PATCH)
	$(WGET) $(CMAKE_TARBALL_URL)
	$(WGET) $(CMAKE_CUDA_PATCH_URL)
	@echo '$(CMAKE_TARBALL_SHA256)  $(CMAKE_TARBALL)' | $(SHA256SUM)
	@echo '$(CMAKE_CUDA_PATCH_SHA256)  $(CMAKE_CUDA_PATCH)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-cmake: .fetch-cmake-$(CMAKE_CUDA_VERSION)

.extract-cmake-$(CMAKE_CUDA_VERSION): .fetch-cmake-$(CMAKE_CUDA_VERSION)
	$(RM) $(CMAKE_BUILD_DIR)
	$(TAR) -xzf $(CMAKE_TARBALL)
	cd $(CMAKE_BUILD_DIR) && $(PATCH) -p1 < $(CURDIR)/$(CMAKE_CUDA_PATCH)
	@$(TOUCH) $@

extract-cmake: .extract-cmake-$(CMAKE_CUDA_VERSION)

.configure-cmake-$(CMAKE_CUDA_VERSION): .extract-cmake-$(CMAKE_CUDA_VERSION)
	cd $(CMAKE_BUILD_DIR) && ./configure --prefix=$(CMAKE_INSTALL_DIR)
	@$(TOUCH) $@

configure-cmake: .configure-cmake-$(CMAKE_CUDA_VERSION)

.build-cmake-$(CMAKE_CUDA_VERSION): .configure-cmake-$(CMAKE_CUDA_VERSION)
	cd $(CMAKE_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-cmake: .build-cmake-$(CMAKE_CUDA_VERSION)

install-cmake: .build-cmake-$(CMAKE_CUDA_VERSION)
	cd $(CMAKE_BUILD_DIR) && make install

clean-cmake:
	@$(RM) .build-cmake-$(CMAKE_CUDA_VERSION)
	@$(RM) .configure-cmake-$(CMAKE_CUDA_VERSION)
	@$(RM) .extract-cmake-$(CMAKE_CUDA_VERSION)
	$(RM) $(CMAKE_BUILD_DIR)

distclean-cmake: clean-cmake
	@$(RM) .fetch-cmake-$(CMAKE_CUDA_VERSION)
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
LUA_VERSION = 5.2.1
LUA_TARBALL_SHA256 = 64304da87976133196f9e4c15250b70f444467b6ed80d7cfd7b3b982b5177be5
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
	cd $(LUA_BUILD_DIR) && make linux CFLAGS="$(LUA_CFLAGS)" LIBS="$(LUA_LIBS)"
	@$(TOUCH) $@

build-lua: .build-lua-$(LUA_VERSION)

install-lua: .build-lua-$(LUA_VERSION)
	cd $(LUA_BUILD_DIR) && make install INSTALL_TOP=$(LUA_INSTALL_DIR)

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

LUAJIT_VERSION = 2.0.0-beta10
LUAJIT_TARBALL = LuaJIT-$(LUAJIT_VERSION).tar.gz
LUAJIT_TARBALL_URL = http://luajit.org/download/$(LUAJIT_TARBALL)
LUAJIT_TARBALL_SHA256 = 750e9317ca2c61fa17e739abc582c55c1fe69749ba65442dfd7f04ce20cf9ff6
LUAJIT_BUILD_DIR = LuaJIT-$(LUAJIT_VERSION)
LUAJIT_INSTALL_DIR = $(PREFIX)/luajit-$(LUAJIT_VERSION)
LUAJIT_CFLAGS = -fPIC

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
	cd $(LUAJIT_BUILD_DIR) && make amalg CFLAGS=$(LUAJIT_CFLAGS) $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-luajit: .build-luajit-$(LUAJIT_VERSION)

install-luajit: .build-luajit-$(LUAJIT_VERSION)
	cd $(LUAJIT_BUILD_DIR) && make install PREFIX=$(LUAJIT_INSTALL_DIR)
	ln -sf luajit-$(LUAJIT_VERSION) $(LUAJIT_INSTALL_DIR)/bin/lua
	ln -sf libluajit-5.1.a $(LUAJIT_INSTALL_DIR)/lib/liblua.a
	ln -sf luajit-2.0 $(LUAJIT_INSTALL_DIR)/include/lua

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
## Boost C++ libraries with Boost.Log
##

BOOST_VERSION = 1.51.0
BOOST_RELEASE = 1_51_0
BOOST_TARBALL = boost_$(BOOST_RELEASE).tar.bz2
BOOST_TARBALL_URL = http://sourceforge.net/projects/boost/files/boost/$(BOOST_VERSION)/$(BOOST_TARBALL)
BOOST_TARBALL_SHA256 = fb2d2335a29ee7fe040a197292bfce982af84a645c81688a915c84c925b69696
BOOST_BUILD_DIR = boost_$(BOOST_RELEASE)
BOOST_INSTALL_DIR = $(PREFIX)/boost_$(BOOST_RELEASE)
BOOST_BUILD_FLAGS = threading=multi variant=release --layout=tagged toolset=gcc cxxflags="-fPIC -std=c++11" dll-path=$(BOOST_INSTALL_DIR)/lib

ifndef USE_BZIP2
BOOST_BUILD_FLAGS += -sNO_BZIP2=1
endif
ifndef USE_PYTHON
BOOST_BUILD_FLAGS += --without-python
endif

define BOOST_PATCH
--- boost/random/detail/integer_log2.hpp
+++ boost/random/detail/integer_log2.hpp
@@ -27,7 +27,7 @@
 #elif defined(BOOST_MSVC)
 #define BOOST_RANDOM_DETAIL_CONSTEXPR __forceinline
 #elif defined(__GNUC__) && __GNUC__ >= 4
-#define BOOST_RANDOM_DETAIL_CONSTEXPR __attribute__((const)) __attribute__((always_inline))
+#define BOOST_RANDOM_DETAIL_CONSTEXPR inline __attribute__((const)) __attribute__((always_inline))
 #else
 #define BOOST_RANDOM_DETAIL_CONSTEXPR inline
 #endif
--- boost/test/detail/global_typedef.hpp
+++ boost/test/detail/global_typedef.hpp
@@ -67,12 +67,13 @@
 // helper templates to prevent ODR violations 
 template<class T> 
 struct static_constant { 
-    static T value; 
+    static T const& value()
+    {
+        static T const v = {};
+        return v;
+    }
 }; 
 
-template<class T> 
-T static_constant<T>::value; 
-
 //____________________________________________________________________________// 
 
 } // namespace ut_detail
--- boost/test/floating_point_comparison.hpp
+++ boost/test/floating_point_comparison.hpp
@@ -248,7 +248,7 @@
 };
 
 namespace {
-check_is_close_t const& check_is_close = unit_test::ut_detail::static_constant<check_is_close_t>::value;
+check_is_close_t const& check_is_close = unit_test::ut_detail::static_constant<check_is_close_t>::value();
 }
 
 //____________________________________________________________________________//
@@ -270,7 +270,7 @@
 };
 
 namespace {
-check_is_small_t const& check_is_small = unit_test::ut_detail::static_constant<check_is_small_t>::value;
+check_is_small_t const& check_is_small = unit_test::ut_detail::static_constant<check_is_small_t>::value();
 }
 
 //____________________________________________________________________________//
--- boost/parameter/keyword.hpp
+++ boost/parameter/keyword.hpp
@@ -91,18 +91,13 @@
     // every instantiation of a function template is the same object.
     // We provide a reference to a common instance of each keyword
     // object and prevent construction by users.
-    static keyword<Tag> const instance;
-
-    // This interface is deprecated
-    static keyword<Tag>& get()
+    static keyword<Tag> const& get()
     {
-        return const_cast<keyword<Tag>&>(instance);
+        static keyword<Tag> const instance = {};
+        return instance;
     }
 };
 
-template <class Tag>
-keyword<Tag> const keyword<Tag>::instance = {};
-
 // Reduces boilerplate required to declare and initialize keywords
 // without violating ODR.  Declares a keyword tag type with the given
 // name in namespace tag_namespace, and declares and initializes a
@@ -123,7 +118,7 @@
       };                                                                \ 
     }                                                                   \ 
     static ::boost::parameter::keyword<tag_namespace::name> const& name \ 
-       = ::boost::parameter::keyword<tag_namespace::name>::instance;
+       = ::boost::parameter::keyword<tag_namespace::name>::get();
 
 #else
 
@@ -141,7 +136,7 @@
     namespace                                                       \ 
     {                                                               \ 
        ::boost::parameter::keyword<tag_namespace::name> const& name \ 
-       = ::boost::parameter::keyword<tag_namespace::name>::instance;\ 
+       = ::boost::parameter::keyword<tag_namespace::name>::get();   \ 
     }
 
 #endif
--- boost/parameter/name.hpp
+++ boost/parameter/name.hpp
@@ -78,13 +78,13 @@
 # if BOOST_WORKAROUND(BOOST_MSVC, < 1300)
 #  define BOOST_PARAMETER_NAME_OBJECT(tag, name)                    \ 
     static ::boost::parameter::keyword<tag> const& name             \ 
-       = ::boost::parameter::keyword<tag>::instance;
+       = ::boost::parameter::keyword<tag>::get();
 # else
 #  define BOOST_PARAMETER_NAME_OBJECT(tag, name)                    \ 
     namespace                                                       \ 
     {                                                               \ 
        ::boost::parameter::keyword<tag> const& name                 \ 
-       = ::boost::parameter::keyword<tag>::instance;                \ 
+       = ::boost::parameter::keyword<tag>::get();                   \ 
     }
 # endif
 
--- boost/parameter/preprocessor.hpp
+++ boost/parameter/preprocessor.hpp
@@ -691,7 +691,7 @@
 # define BOOST_PARAMETER_FUNCTION_DEFAULT_EVAL_DEFAULT(arg, tag_namespace) \ 
     boost::parameter::keyword< \ 
         tag_namespace::BOOST_PARAMETER_FN_ARG_KEYWORD(arg) \ 
-    >::instance | boost::parameter::aux::use_default_tag()
+    >::get() | boost::parameter::aux::use_default_tag()
 
 # define BOOST_PARAMETER_FUNCTION_DEFAULT_FUNCTION_GET_ARG(arg, tag_ns) \ 
     BOOST_PARAMETER_FUNCTION_CAST( \ 
@@ -842,7 +842,7 @@
 # define BOOST_PARAMETER_FUNCTION_DEFAULT_GET_ARG(r, tag_ns, arg) \ 
     , BOOST_PARAMETER_FUNCTION_CAST( \ 
           args[ \ 
-              boost::parameter::keyword<tag_ns::BOOST_PARAMETER_FN_ARG_KEYWORD(arg)>::instance \ 
+              boost::parameter::keyword<tag_ns::BOOST_PARAMETER_FN_ARG_KEYWORD(arg)>::get() \ 
           ] \ 
         , BOOST_PARAMETER_FN_ARG_PRED(arg) \ 
         , Args \ 
endef
export BOOST_PATCH

.fetch-boost-$(BOOST_VERSION):
	@$(RM) $(BOOST_TARBALL)
	$(WGET) $(BOOST_TARBALL_URL)
	@echo '$(BOOST_TARBALL_SHA256)  $(BOOST_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-boost: .fetch-boost-$(BOOST_VERSION)

.extract-boost-$(BOOST_VERSION): .fetch-boost-$(BOOST_VERSION)
	$(TAR) -xjf $(BOOST_TARBALL)
	cd $(BOOST_BUILD_DIR) && echo "$$BOOST_PATCH" | sed -e 's/\\ $$/\\/' | $(PATCH) -p0
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
	@echo
	@echo '# add Boost $(BOOST_VERSION) to environment'
	@echo 'export PATH="$(BOOST_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export LD_LIBRARY_PATH="$(BOOST_INSTALL_DIR)/lib$${LD_LIBRARY_PATH+:$$LD_LIBRARY_PATH}"'
	@echo 'export PYTHONPATH="$(BOOST_INSTALL_DIR)/lib$${PYTHONPATH+:$$PYTHONPATH}"'
	@echo 'export CMAKE_PREFIX_PATH="$(BOOST_INSTALL_DIR)$${CMAKE_PREFIX_PATH+:$$CMAKE_PREFIX_PATH}"'

##
## HDF5 C++ library
##

HDF5_VERSION = 1.8.9
HDF5_TARBALL = hdf5-$(HDF5_VERSION).tar.bz2
HDF5_TARBALL_URL = http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_VERSION)/src/$(HDF5_TARBALL)
HDF5_TARBALL_SHA256 = fafe54856b00f0d6531629f66e2c476ab6ee03458803088f270bbfc4a60966c7
HDF5_BUILD_DIR = hdf5-$(HDF5_VERSION)
HDF5_INSTALL_DIR = $(PREFIX)/hdf5-$(HDF5_VERSION)
HDF5_CONFIGURE_FLAGS = --enable-cxx
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
	cd $(HDF5_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-hdf5: .build-hdf5-$(HDF5_VERSION)

install-hdf5: .build-hdf5-$(HDF5_VERSION)
	cd $(HDF5_BUILD_DIR) && make install

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

GIT_VERSION = 1.7.11.4
GIT_TARBALL = git-$(GIT_VERSION).tar.gz
GIT_TARBALL_URL = http://git-core.googlecode.com/files/$(GIT_TARBALL)
GIT_TARBALL_SHA256 = e91d6d2552d2def0e1431feb2e521fb40c0722f6ceff965592027153447d4199
GIT_MANPAGES_TARBALL = git-manpages-$(GIT_VERSION).tar.gz
GIT_MANPAGES_TARBALL_URL = http://git-core.googlecode.com/files/$(GIT_MANPAGES_TARBALL)
GIT_MANPAGES_TARBALL_SHA256 = 0ed00facfb6d1744f31bdb493f1626c1901c930f0eec7cddea93262c7e022b33
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
	cd $(GIT_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-git: .build-git-$(GIT_VERSION)

install-git: .build-git-$(GIT_VERSION)
	cd $(GIT_BUILD_DIR) && make install
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
## htop
##

HTOP_VERSION = 1.0.1
HTOP_TARBALL = htop-$(HTOP_VERSION).tar.gz
HTOP_TARBALL_URL = http://sourceforge.net/projects/htop/files/htop/$(HTOP_VERSION)/$(HTOP_TARBALL)
HTOP_TARBALL_SHA256 = 07db2cbe02835f9e186b9610ecc3beca330a5c9beadb3b6069dd0a10561506f2
HTOP_BUILD_DIR = htop-$(HTOP_VERSION)
HTOP_INSTALL_DIR = $(PREFIX)/htop-$(HTOP_VERSION)

.fetch-htop-$(HTOP_VERSION):
	@$(RM) $(HTOP_TARBALL)
	$(WGET) $(HTOP_TARBALL_URL)
	@echo '$(HTOP_TARBALL_SHA256)  $(HTOP_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-htop: .fetch-htop-$(HTOP_VERSION)

.extract-htop-$(HTOP_VERSION): .fetch-htop-$(HTOP_VERSION)
	$(RM) $(HTOP_BUILD_DIR)
	$(TAR) -xzf $(HTOP_TARBALL)
	@$(TOUCH) $@

extract-htop: .extract-htop-$(HTOP_VERSION)

.configure-htop-$(HTOP_VERSION): .extract-htop-$(HTOP_VERSION)
	cd $(HTOP_BUILD_DIR) && ./configure --prefix=$(HTOP_INSTALL_DIR)
	@$(TOUCH) $@

configure-htop: .configure-htop-$(HTOP_VERSION)

.build-htop-$(HTOP_VERSION): .configure-htop-$(HTOP_VERSION)
	cd $(HTOP_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-htop: .build-htop-$(HTOP_VERSION)

install-htop: .build-htop-$(HTOP_VERSION)
	cd $(HTOP_BUILD_DIR) && make install

clean-htop:
	@$(RM) .build-htop-$(HTOP_VERSION)
	@$(RM) .configure-htop-$(HTOP_VERSION)
	@$(RM) .extract-htop-$(HTOP_VERSION)
	$(RM) $(HTOP_BUILD_DIR)

distclean-htop: clean-htop
	@$(RM) .fetch-htop-$(HTOP_VERSION)
	$(RM) $(HTOP_TARBALL)

env-htop:
	@echo
	@echo '# add htop $(HTOP_VERSION) to environment'
	@echo 'export PATH="$(HTOP_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(HTOP_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

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
	cd $(GRAPHVIZ_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-graphviz: .build-graphviz-$(GRAPHVIZ_VERSION)

install-graphviz: .build-graphviz-$(GRAPHVIZ_VERSION)
	cd $(GRAPHVIZ_BUILD_DIR) && make install

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

CLANG_VERSION = 3.1
LLVM_TARBALL = llvm-$(CLANG_VERSION).src.tar.gz
LLVM_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(LLVM_TARBALL)
LLVM_TARBALL_SHA256 = 1ea05135197b5400c1f88d00ff280d775ce778f8f9ea042e25a1e1e734a4b9ab
CLANG_TARBALL = clang-$(CLANG_VERSION).src.tar.gz
CLANG_TARBALL_URL = http://llvm.org/releases/$(CLANG_VERSION)/$(CLANG_TARBALL)
CLANG_TARBALL_SHA256 = ff63e215dcd3e2838ffdea38502f8d35bab17e487f3c3799579961e452d5a786
CLANG_BUILD_DIR = llvm-$(CLANG_VERSION).src
CLANG_CONFIGURE_FLAGS = --enable-optimized --enable-bindings=none --with-gcc-toolchain=$(GCC_INSTALL_DIR)
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
	cd $(CLANG_BUILD_DIR)/tools && $(TAR) -xzf $(CURDIR)/$(CLANG_TARBALL) && mv clang-$(CLANG_VERSION).src clang
	@$(TOUCH) $@

extract-clang: .extract-clang-$(CLANG_VERSION)

.configure-clang-$(CLANG_VERSION): .extract-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && ./configure $(CLANG_CONFIGURE_FLAGS) --prefix=$(CLANG_INSTALL_DIR)
	@$(TOUCH) $@

configure-clang: .configure-clang-$(CLANG_VERSION)

.build-clang-$(CLANG_VERSION): .configure-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-clang: .build-clang-$(CLANG_VERSION)

install-clang: .build-clang-$(CLANG_VERSION)
	cd $(CLANG_BUILD_DIR) && make install

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
	cd $(GNU_PARALLEL_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gnu-parallel: .build-gnu-parallel-$(GNU_PARALLEL_VERSION)

install-gnu-parallel: .build-gnu-parallel-$(GNU_PARALLEL_VERSION)
	cd $(GNU_PARALLEL_BUILD_DIR) && make install

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
GMP_TARBALL_URL = ftp://ftp.gmplib.org/pub/gmp-$(GMP_VERSION)/$(GMP_TARBALL)
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
	cd $(GMP_BUILD_DIR) && CFLAGS=-fPIC ./configure --prefix=$(GMP_INSTALL_DIR) --enable-cxx --disable-shared
	@$(TOUCH) $@

configure-gmp: .configure-gmp-$(GMP_VERSION)

.build-gmp-$(GMP_VERSION): .configure-gmp-$(GMP_VERSION)
	cd $(GMP_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gmp: .build-gmp-$(GMP_VERSION)

.install-gmp-$(GMP_VERSION): .build-gmp-$(GMP_VERSION)
	cd $(GMP_BUILD_DIR) && make install
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
	cd $(MPFR_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-mpfr: .build-mpfr-$(MPFR_VERSION)

.install-mpfr-$(MPFR_VERSION): .build-mpfr-$(MPFR_VERSION)
	cd $(MPFR_BUILD_DIR) && make install
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
	cd $(MPC_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-mpc: .build-mpc-$(MPC_VERSION)

.install-mpc-$(MPC_VERSION): .build-mpc-$(MPC_VERSION)
	cd $(MPC_BUILD_DIR) && make install
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
	cd $(PPL_BUILD_DIR) && CPPFLAGS=-I$(GMP_INSTALL_DIR)/include LDFLAGS=-L$(GMP_INSTALL_DIR)/lib ./configure --prefix=$(PPL_INSTALL_DIR) --enable-interfaces=c,cxx --disable-shared --disable-watchdog
	@$(TOUCH) $@

configure-ppl: .configure-ppl-$(PPL_VERSION)

.build-ppl-$(PPL_VERSION): .configure-ppl-$(PPL_VERSION)
	cd $(PPL_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-ppl: .build-ppl-$(PPL_VERSION)

.install-ppl-$(PPL_VERSION): .build-ppl-$(PPL_VERSION)
	cd $(PPL_BUILD_DIR) && make install
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
	cd $(CLOOG_PPL_BUILD_DIR) && ./configure --prefix=$(CLOOG_PPL_INSTALL_DIR) --with-ppl=$(PPL_INSTALL_DIR) --with-gmp=$(GMP_INSTALL_DIR) --disable-shared --with-host-libstdcxx=-lstdc++
	@$(TOUCH) $@

configure-cloog-ppl: .configure-cloog-ppl-$(CLOOG_PPL_VERSION)

.build-cloog-ppl-$(CLOOG_PPL_VERSION): .configure-cloog-ppl-$(CLOOG_PPL_VERSION)
	cd $(CLOOG_PPL_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-cloog-ppl: .build-cloog-ppl-$(CLOOG_PPL_VERSION)

.install-cloog-ppl-$(CLOOG_PPL_VERSION): .build-cloog-ppl-$(CLOOG_PPL_VERSION)
	cd $(CLOOG_PPL_BUILD_DIR) && make install
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

GCC_VERSION = 4.7.1
GCC_TARBALL = gcc-$(GCC_VERSION).tar.bz2
GCC_TARBALL_URL = http://ftp.gwdg.de/pub/misc/gcc/releases/gcc-$(GCC_VERSION)/$(GCC_TARBALL)
GCC_TARBALL_SHA256 = 16093f6fa01732adf378d97fe338f113c933bdf56da22bf87c76beff13da406f
GCC_BUILD_DIR = gcc-$(GCC_VERSION)
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
	cd $(GCC_BUILD_DIR) && LDFLAGS=-L$(GMP_INSTALL_DIR)/lib ac_cv_lib_pwl_PWL_handle_timeout=no ./configure --prefix=$(GCC_INSTALL_DIR) --disable-multilib --enable-languages=c,c++,fortran,lto --with-gmp=$(GMP_INSTALL_DIR) --with-mpfr=$(MPFR_INSTALL_DIR) --with-mpc=$(MPC_INSTALL_DIR) --with-ppl=$(PPL_INSTALL_DIR) --with-cloog=$(CLOOG_PPL_INSTALL_DIR) --enable-build-with-cxx
	@$(TOUCH) $@

configure-gcc: .configure-gcc-$(GCC_VERSION)

.build-gcc-$(GCC_VERSION): .configure-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-gcc: .build-gcc-$(GCC_VERSION)

install-gcc: .build-gcc-$(GCC_VERSION)
	cd $(GCC_BUILD_DIR) && make install

clean-gcc:
	@$(RM) .build-gcc-$(GCC_VERSION)
	@$(RM) .configure-gcc-$(GCC_VERSION)
	@$(RM) .extract-gcc-$(GCC_VERSION)
	@$(RM) .install-gmp-$(GMP_VERSION)
	@$(RM) .install-mpfr-$(MPFR_VERSION)
	@$(RM) .install-mpc-$(MPC_VERSION)
	@$(RM) .install-ppl-$(PPL_VERSION)
	@$(RM) .install-cloog-ppl-$(CLOOG_PPL_VERSION)
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
## HALMD Highly Accerelated Large-Scale Molecular Dynamics
##

HALMD_VERSION = 0.2-rc4
HALMD_TARBALL = halmd-$(HALMD_VERSION).tar.bz2
HALMD_TARBALL_URL = http://downloads.sourceforge.net/project/halmd/halmd/0.2/testing/$(HALMD_TARBALL)
HALMD_TARBALL_SHA256 = 20e4405c15e4848f76a4afd4b8d94c7aada1eddc05811ec74e6c59d315f77189
HALMD_SOURCE_DIR = halmd-$(HALMD_VERSION)
HALMD_BUILD_DIR = $(HALMD_SOURCE_DIR)/build/release
HALMD_INSTALL_DIR = $(PREFIX)/halmd-$(HALMD_VERSION)
HALMD_BUILD_ENV = CUDACC="nvcc --compiler-bindir=/usr/bin"

.fetch-halmd-$(HALMD_VERSION):
	@$(RM) $(HALMD_TARBALL)
	$(WGET) $(HALMD_TARBALL_URL)
	@echo '$(HALMD_TARBALL_SHA256)  $(HALMD_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-halmd: .fetch-halmd-$(HALMD_VERSION)

.extract-halmd-$(HALMD_VERSION): .fetch-halmd-$(HALMD_VERSION)
	$(RM) $(HALMD_SOURCE_DIR)
	$(TAR) -xjf $(HALMD_TARBALL)
	@$(TOUCH) $@

extract-halmd: .extract-halmd-$(HALMD_VERSION)

.configure-halmd-$(HALMD_VERSION): .extract-halmd-$(HALMD_VERSION)
	mkdir -p $(HALMD_BUILD_DIR)
	cd $(HALMD_BUILD_DIR) && $(HALMD_BUILD_ENV) cmake -DCMAKE_INSTALL_PREFIX=$(HALMD_INSTALL_DIR) $(CURDIR)/$(HALMD_SOURCE_DIR)
	@$(TOUCH) $@

configure-halmd: .configure-halmd-$(HALMD_VERSION)

.build-halmd-$(HALMD_VERSION): .configure-halmd-$(HALMD_VERSION)
	cd $(HALMD_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-halmd: .build-halmd-$(HALMD_VERSION)

install-halmd: .build-halmd-$(HALMD_VERSION)
	cd $(HALMD_BUILD_DIR) && make install

clean-halmd:
	@$(RM) .build-halmd-$(HALMD_VERSION)
	@$(RM) .configure-halmd-$(HALMD_VERSION)
	@$(RM) .extract-halmd-$(HALMD_VERSION)
	$(RM) $(HALMD_SOURCE_DIR)

distclean-halmd: clean-halmd
	@$(RM) .fetch-halmd-$(HALMD_VERSION)
	$(RM) $(HALMD_TARBALL)

env-halmd:
	@echo
	@echo '# add CMake $(HALMD_VERSION) to environment'
	@echo 'export PATH="$(HALMD_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(HALMD_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'

##
## nvCUDA-tools - A collection of tools for NVIDIA CUDA compute devices
##

NVCUDA_TOOLS_VERSION = 15490f0
NVCUDA_TOOLS_GIT_URL = http://git.colberg.org/nvcuda-tools.git
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
	mkdir -p $(NVCUDA_TOOLS_BUILD_DIR)
	cd $(NVCUDA_TOOLS_BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(NVCUDA_TOOLS_INSTALL_DIR) $(CURDIR)/$(NVCUDA_TOOLS_SOURCE_DIR)
	@$(TOUCH) $@

configure-nvcuda-tools: .configure-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)

.build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION): .configure-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	cd $(NVCUDA_TOOLS_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-nvcuda-tools: .build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)

install-nvcuda-tools: .build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	cd $(NVCUDA_TOOLS_BUILD_DIR) && make install

clean-nvcuda-tools:
	@$(RM) .build-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	@$(RM) .configure-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	$(RM) $(NVCUDA_TOOLS_BUILD_DIR)

distclean-nvcuda-tools: clean-nvcuda-tools
	@$(RM) .fetch-nvcuda-tools-$(NVCUDA_TOOLS_VERSION)
	$(RM) $(NVCUDA_TOOLS_SOURCE_DIR)

env-nvcuda-tools:
	@echo
	@echo '# add CMake $(NVCUDA_TOOLS_VERSION) to environment'
	@echo 'export PATH="$(NVCUDA_TOOLS_INSTALL_DIR)/bin$${PATH+:$$PATH}"'

##
## Bootstraping GHC (Glasgow Haskell Compiler)
##

GHC_BOOTSTRAP_VERSION = 6.12.3
GHC_BOOTSTRAP_TARBALL = ghc-$(GHC_BOOTSTRAP_VERSION)-x86_64-unknown-linux-n.tar.bz2
GHC_BOOTSTRAP_TARBALL_URL =  http://www.haskell.org/ghc/dist/$(GHC_BOOTSTRAP_VERSION)/$(GHC_BOOTSTRAP_TARBALL)
GHC_BOOTSTRAP_TARBALL_SHA256 = f0e13bdec040f06b1074595ddc39064f75214aee05d64554f4809dca3e138001
GHC_BOOTSTRAP_BUILD_DIR = ghc-$(GHC_BOOTSTRAP_VERSION)
GHC_BOOTSTRAP_INSTALL_DIR = $(CURDIR)/.ghc-$(GHC_BOOTSTRAP_VERSION)

.fetch-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION):
	@$(RM) $(GHC_BOOTSTRAP_TARBALL)
	$(WGET) $(GHC_BOOTSTRAP_TARBALL_URL)
	@echo '$(GHC_BOOTSTRAP_TARBALL_SHA256)  $(GHC_BOOTSTRAP_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-ghc-bootstrap: .fetch-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)

.extract-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION): .fetch-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)
	$(RM) $(GHC_BOOTSTRAP_BUILD_DIR)
	$(TAR) -xjf $(GHC_BOOTSTRAP_TARBALL)
	@$(TOUCH) $@

extract-ghc-bootstrap: .extract-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)

.configure-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION): .extract-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)
	cd $(GHC_BOOTSTRAP_BUILD_DIR) && ./configure --prefix=$(GHC_BOOTSTRAP_INSTALL_DIR)
	@$(TOUCH) $@

configure-ghc-bootstrap: .configure-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)

.install-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION): .configure-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION) .install-gmp-$(GMP_VERSION)
	cd $(GHC_BOOTSTRAP_BUILD_DIR) && make install
	sed -i -e 's/-pgmc "\([^"]*\)"/-pgmc "\1" -pgma "\1" -pgml "\1" -pgmP "\1 -E -undef -traditional"/' $(GHC_BOOTSTRAP_INSTALL_DIR)/bin/ghc-$(GHC_BOOTSTRAP_VERSION)
	ln -s -f -t $(GHC_BOOTSTRAP_INSTALL_DIR)/lib/ghc-$(GHC_BOOTSTRAP_VERSION) $(GMP_INSTALL_DIR)/lib/libgmp.a
	@$(TOUCH) $@

clean-ghc-bootstrap:
	@$(RM) .configure-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)
	@$(RM) .extract-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)
	$(RM) $(GHC_BOOTSTRAP_BUILD_DIR)

distclean-ghc-bootstrap: clean-ghc-bootstrap
	@$(RM) .fetch-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION)
	$(RM) $(GHC_BOOTSTRAP_TARBALL)

##
## GHC (Glasgow Haskell Compiler)
##

GHC_VERSION = 7.4.2
GHC_TARBALL = ghc-$(GHC_VERSION)-src.tar.bz2
GHC_TARBALL_URL =  http://www.haskell.org/ghc/dist/$(GHC_VERSION)/$(GHC_TARBALL)
GHC_TARBALL_SHA256 = f2ee1289a33cc70539287129841acc7eaf16112bb60c59b5a6ee91887bfd836d
GHC_BUILD_DIR = ghc-$(GHC_VERSION)
GHC_INSTALL_DIR = $(CURDIR)/.ghc-$(GHC_VERSION)

.fetch-ghc-$(GHC_VERSION):
	@$(RM) $(GHC_TARBALL)
	$(WGET) $(GHC_TARBALL_URL)
	@echo '$(GHC_TARBALL_SHA256)  $(GHC_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-ghc: .fetch-ghc-$(GHC_VERSION)

.extract-ghc-$(GHC_VERSION): .fetch-ghc-$(GHC_VERSION)
	$(RM) $(GHC_BUILD_DIR)
	$(TAR) -xjf $(GHC_TARBALL)
	@$(TOUCH) $@

extract-ghc: .extract-ghc-$(GHC_VERSION)

.configure-ghc-$(GHC_VERSION): .extract-ghc-$(GHC_VERSION) .install-ghc-bootstrap-$(GHC_BOOTSTRAP_VERSION) .install-gmp-$(GMP_VERSION)
	cd $(GHC_BUILD_DIR) && LIBRARY_PATH=$(GMP_INSTALL_DIR)/lib PATH=$(GHC_BOOTSTRAP_INSTALL_DIR)/bin:$(PATH) ./configure --prefix=$(GHC_INSTALL_DIR) --with-gmp-includes=$(GMP_INSTALL_DIR)/include --with-gmp-libraries=$(GMP_INSTALL_DIR)/lib
	@$(TOUCH) $@

configure-ghc: .configure-ghc-$(GHC_VERSION)

.build-ghc-$(GHC_VERSION): .configure-ghc-$(GHC_VERSION)
	cd $(GHC_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-ghc: .build-ghc-$(GHC_VERSION)

.install-ghc-$(GHC_VERSION): .build-ghc-$(GHC_VERSION)
	cd $(GHC_BUILD_DIR) && make install
	@$(TOUCH) $@

clean-ghc:
	@$(RM) .configure-ghc-$(GHC_VERSION)
	@$(RM) .build-ghc-$(GHC_VERSION)
	@$(RM) .extract-ghc-$(GHC_VERSION)
	$(RM) $(GHC_BUILD_DIR)

distclean-ghc: clean-ghc
	@$(RM) .fetch-ghc-$(GHC_VERSION)
	$(RM) $(GHC_TARBALL)

##
## Haskell Cabal
##

HASKELL_CABAL_VERSION = 0.14.0
HASKELL_CABAL_TARBALL = cabal-install-$(HASKELL_CABAL_VERSION).tar.gz
HASKELL_CABAL_TARBALL_URL = http://hackage.haskell.org/packages/archive/cabal-install/$(HASKELL_CABAL_VERSION)/$(HASKELL_CABAL_TARBALL)
HASKELL_CABAL_TARBALL_SHA256 = f4f2b50269ff59d67b5f3d82d50f7706b6bad1117295a7c81f32bbe72add5bd8
HASKELL_CABAL_BUILD_DIR = cabal-install-$(HASKELL_CABAL_VERSION)
HASKELL_CABAL_INSTALL_DIR = $(CURDIR)/.cabal-install-$(HASKELL_CABAL_VERSION)

.fetch-haskell-cabal-$(HASKELL_CABAL_VERSION):
	@$(RM) $(HASKELL_CABAL_TARBALL)
	$(WGET) $(HASKELL_CABAL_TARBALL_URL)
	@echo '$(HASKELL_CABAL_TARBALL_SHA256)  $(HASKELL_CABAL_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-haskell-cabal: .fetch-haskell-cabal-$(HASKELL_CABAL_VERSION)

.extract-haskell-cabal-$(HASKELL_CABAL_VERSION): .fetch-haskell-cabal-$(HASKELL_CABAL_VERSION)
	$(RM) $(HASKELL_CABAL_BUILD_DIR)
	$(TAR) -xzf $(HASKELL_CABAL_TARBALL)
	@$(TOUCH) $@

extract-haskell-cabal: .extract-haskell-cabal-$(HASKELL_CABAL_VERSION)

.install-haskell-cabal-$(HASKELL_CABAL_VERSION): .extract-haskell-cabal-$(HASKELL_CABAL_VERSION) .install-ghc-$(GHC_VERSION)
	cd $(HASKELL_CABAL_BUILD_DIR) && HOME=$(CURDIR)/$(HASKELL_CABAL_BUILD_DIR) PREFIX=$(HASKELL_CABAL_INSTALL_DIR) PATH=$(GHC_INSTALL_DIR)/bin:$(PATH) sh bootstrap.sh
	@$(TOUCH) $@

clean-haskell-cabal:
	@$(RM) .extract-haskell-cabal-$(HASKELL_CABAL_VERSION)
	$(RM) $(HASKELL_CABAL_BUILD_DIR)

distclean-haskell-cabal: clean-haskell-cabal
	@$(RM) .fetch-haskell-cabal-$(HASKELL_CABAL_VERSION)
	$(RM) $(HASKELL_CABAL_TARBALL)

##
## PCRE (Perl Compatible Regular Expressions)
##

PCRE_VERSION = 8.30
PCRE_TARBALL = pcre-$(PCRE_VERSION).tar.bz2
PCRE_TARBALL_URL = http://sourceforge.net/projects/pcre/files/pcre/$(PCRE_VERSION)/$(PCRE_TARBALL)
PCRE_TARBALL_SHA256 = c1113fd7db934e97ad8b3917d432e5b642e9eb9afd127eb797804937c965f4ac
PCRE_BUILD_DIR = pcre-$(PCRE_VERSION)
PCRE_INSTALL_DIR = $(CURDIR)/.pcre-$(PCRE_VERSION)

.fetch-pcre-$(PCRE_VERSION):
	@$(RM) $(PCRE_TARBALL)
	$(WGET) $(PCRE_TARBALL_URL)
	@echo '$(PCRE_TARBALL_SHA256)  $(PCRE_TARBALL)' | $(SHA256SUM)
	@$(TOUCH) $@

fetch-pcre: .fetch-pcre-$(PCRE_VERSION)

.extract-pcre-$(PCRE_VERSION): .fetch-pcre-$(PCRE_VERSION)
	$(RM) $(PCRE_BUILD_DIR)
	$(TAR) -xjf $(PCRE_TARBALL)
	@$(TOUCH) $@

extract-pcre: .extract-pcre-$(PCRE_VERSION)

.configure-pcre-$(PCRE_VERSION): .extract-pcre-$(PCRE_VERSION)
	cd $(PCRE_BUILD_DIR) && CFLAGS=-fPIC ./configure --prefix=$(PCRE_INSTALL_DIR) --disable-shared
	@$(TOUCH) $@

configure-pcre: .configure-pcre-$(PCRE_VERSION)

.build-pcre-$(PCRE_VERSION): .configure-pcre-$(PCRE_VERSION)
	cd $(PCRE_BUILD_DIR) && make $(PARALLEL_BUILD_FLAGS)
	@$(TOUCH) $@

build-pcre: .build-pcre-$(PCRE_VERSION)

.install-pcre-$(PCRE_VERSION): .build-pcre-$(PCRE_VERSION)
	cd $(PCRE_BUILD_DIR) && make install
	@$(TOUCH) $@

clean-pcre:
	@$(RM) .build-pcre-$(PCRE_VERSION)
	@$(RM) .configure-pcre-$(PCRE_VERSION)
	@$(RM) .extract-pcre-$(PCRE_VERSION)
	$(RM) $(PCRE_BUILD_DIR)

distclean-pcre: clean-pcre
	@$(RM) .fetch-pcre-$(PCRE_VERSION)
	$(RM) $(PCRE_TARBALL)

##
## git-annex
##

GIT_ANNEX_VERSION = 3.20120807
GIT_ANNEX_BUILD_DIR = git-annex-$(GIT_ANNEX_VERSION)
GIT_ANNEX_INSTALL_DIR = $(PREFIX)/git-annex-$(GIT_ANNEX_VERSION)

.configure-git-annex-$(GIT_ANNEX_VERSION): .install-ghc-$(GHC_VERSION) .install-haskell-cabal-$(HASKELL_CABAL_VERSION)
	HOME=$(CURDIR)/$(GIT_ANNEX_BUILD_DIR) PATH=$(GHC_INSTALL_DIR)/bin:$(HASKELL_CABAL_INSTALL_DIR)/bin:$(PATH) cabal update
	@$(TOUCH) $@

configure-git-annex: .configure-git-annex-$(GIT_ANNEX_VERSION)

.build-git-annex-$(GIT_ANNEX_VERSION): .configure-git-annex-$(GIT_ANNEX_VERSION) .install-pcre-$(PCRE_VERSION)
	HOME=$(CURDIR)/$(GIT_ANNEX_BUILD_DIR) PATH=$(GHC_INSTALL_DIR)/bin:$(HASKELL_CABAL_INSTALL_DIR)/bin:$(PATH) cabal install git-annex-$(GIT_ANNEX_VERSION) --extra-include-dirs=$(PCRE_INSTALL_DIR)/include --extra-lib-dirs=$(PCRE_INSTALL_DIR)/lib --only-dependencies
	@$(TOUCH) $@

build-git-annex: .build-git-annex-$(GIT_ANNEX_VERSION)

install-git-annex: .build-git-annex-$(GIT_ANNEX_VERSION)
	HOME=$(CURDIR)/$(GIT_ANNEX_BUILD_DIR) PATH=$(GHC_INSTALL_DIR)/bin:$(HASKELL_CABAL_INSTALL_DIR)/bin:$(PATH) cabal install git-annex-$(GIT_ANNEX_VERSION) --extra-include-dirs=$(PCRE_INSTALL_DIR)/include --extra-lib-dirs=$(PCRE_INSTALL_DIR)/lib --prefix=$(GIT_ANNEX_INSTALL_DIR)

clean-git-annex:
	@$(RM) .configure-git-annex-$(GIT_ANNEX_VERSION)
	@$(RM) .build-git-annex-$(GIT_ANNEX_VERSION)
	$(RM) $(GIT_ANNEX_BUILD_DIR)

env-git-annex:
	@echo
	@echo '# add git-annex $(GIT_ANNEX_VERSION) to environment'
	@echo 'export PATH="$(GIT_ANNEX_INSTALL_DIR)/bin$${PATH+:$$PATH}"'
	@echo 'export MANPATH="$(GIT_ANNEX_INSTALL_DIR)/share/man$${MANPATH+:$$MANPATH}"'
