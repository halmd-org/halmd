/* dynamic linking loader wrapper
 *
 * Copyright © 2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <boost/filesystem.hpp>
#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <dlfcn.h>
#include <string>

namespace dl
{

class error : public std::runtime_error
{
public:
    error(char const* str) : std::runtime_error(str) {}
};

class library : boost::noncopyable
{
public:
    library() : handle(NULL) {}
    library(boost::filesystem::path path) { open(path); }
    ~library() { close(); }

    void open(boost::filesystem::path path)
    {
        if (handle)
            close();

        if (path.has_parent_path() && !boost::filesystem::exists(path)) {
            // fall back to default search paths, as described in dlopen(3)
            path = path.filename();
        }
        handle = dlopen(path.string().c_str(), RTLD_GLOBAL|RTLD_NOW);
        if (!handle) {
            throw error(dlerror());
        }
    }

    void close()
    {
        if (handle)
            dlclose(handle);
        handle = NULL;
    }

    template <typename T>
    T* symbol(std::string const& name)
    {
        // reset errors
        dlerror();
        T* sym = reinterpret_cast<T*>(dlsym(handle, name.c_str()));
        char const* dlsym_error = dlerror();
        if (dlsym_error)
            throw error(dlsym_error);
        return sym;
    }

protected:
    void* handle;
};

template <typename T>
class symbol : public boost::function<T>
{
public:
    typedef boost::function<T> _Base;

public:
    symbol(library& dl, std::string const& name) : _Base(dl.symbol<T>(name)) {}
    symbol() {}

    void set(library& dl, std::string const& name)
    {
        _Base::operator=(dl.symbol<T>(name));
    }
};

} // namespace dl
