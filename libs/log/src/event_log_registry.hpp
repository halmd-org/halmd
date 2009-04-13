/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 * 
 * \file   event_log_registry.hpp
 * \author Andrey Semashev
 * \date   16.11.2008
 * 
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#ifndef BOOST_LOG_EVENT_LOG_REGISTRY_HPP_INCLUDED_
#define BOOST_LOG_EVENT_LOG_REGISTRY_HPP_INCLUDED_

#include <windows.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <boost/version.hpp>
#include <boost/optional.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/detail/code_conversion.hpp>
#include <boost/log/detail/throw_exception.hpp>

#ifdef _MSC_VER
#endif

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sinks {

namespace aux {

    //! Helper traits to integrate with WinAPI
    template< typename CharT >
    struct registry_traits;

#ifdef BOOST_LOG_USE_CHAR
    template< >
    struct registry_traits< char >
    {
        static std::string make_event_log_key(std::string const& log_name, std::string const& source_name)
        {
            return "SYSTEM\\CurrentControlSet\\Services\\EventLog\\" + log_name + "\\" + source_name;
        }

        static std::string make_default_log_name()
        {
            return "Application";
        }

        static std::string make_default_source_name()
        {
            char buf[MAX_PATH];
            DWORD size = GetModuleFileNameA(NULL, buf, sizeof(buf) / sizeof(*buf));

            std::string source_name(buf, buf + size);
            if (source_name.empty())
            {
                // In case of error we provide artifical application name
                std::ostringstream strm;
                strm << "Boost.Log "
                    << static_cast< unsigned int >(BOOST_VERSION / 100000)
                    << "."
                    << static_cast< unsigned int >(BOOST_VERSION / 100 % 1000)
                    << "."
                    << static_cast< unsigned int >(BOOST_VERSION % 100);
                source_name = strm.str();
            }
            else
            {
                // Cut off the path and extension
                std::size_t backslash_pos = source_name.rfind('\\');
                if (backslash_pos == std::string::npos || backslash_pos >= source_name.size() - 1)
                    backslash_pos = 0;
                else
                    ++backslash_pos;
                std::size_t dot_pos = source_name.rfind('.');
                if (dot_pos == std::string::npos || dot_pos < backslash_pos)
                    dot_pos = source_name.size();
                source_name = source_name.substr(backslash_pos, dot_pos - backslash_pos);
            }

            return source_name;
        }

        static LSTATUS create_key(
            HKEY hKey,
            const char* lpSubKey,
            DWORD Reserved,
            char* lpClass,
            DWORD dwOptions,
            REGSAM samDesired,
            LPSECURITY_ATTRIBUTES lpSecurityAttributes,
            PHKEY phkResult,
            LPDWORD lpdwDisposition)
        {
            return RegCreateKeyExA(hKey, lpSubKey, Reserved, lpClass, dwOptions, samDesired, lpSecurityAttributes, phkResult, lpdwDisposition);
        }

        static LSTATUS set_value(
            HKEY hKey,
            const char* lpValueName,
            DWORD Reserved,
            DWORD dwType,
            const BYTE *lpData,
            DWORD cbData)
        {
            return RegSetValueExA(hKey, lpValueName, Reserved, dwType, lpData, cbData);
        }

        static const char* get_event_message_file_param_name() { return "EventMessageFile"; }
        static const char* get_category_message_file_param_name() { return "CategoryMessageFile"; }
        static const char* get_category_count_param_name() { return "CategoryCount"; }
        static const char* get_types_supported_param_name() { return "TypesSupported"; }
    };
#endif // BOOST_LOG_USE_CHAR

#ifdef BOOST_LOG_USE_WCHAR_T
    template< >
    struct registry_traits< wchar_t >
    {
        static std::wstring make_event_log_key(std::wstring const& log_name, std::wstring const& source_name)
        {
            return L"SYSTEM\\CurrentControlSet\\Services\\EventLog\\" + log_name + L"\\" + source_name;
        }

        static std::wstring make_default_log_name()
        {
            return L"Application";
        }

        static std::wstring make_default_source_name()
        {
            wchar_t buf[MAX_PATH];
            DWORD size = GetModuleFileNameW(NULL, buf, sizeof(buf) / sizeof(*buf));

            std::wstring source_name(buf, buf + size);
            if (source_name.empty())
            {
                // In case of error we provide artifical application name
                std::wostringstream strm;
                strm << L"Boost.Log "
                    << static_cast< unsigned int >(BOOST_VERSION / 100000)
                    << L"."
                    << static_cast< unsigned int >(BOOST_VERSION / 100 % 1000)
                    << L"."
                    << static_cast< unsigned int >(BOOST_VERSION % 100);
                source_name = strm.str();
            }
            else
            {
                // Cut off the path and extension
                std::size_t backslash_pos = source_name.rfind(L'\\');
                if (backslash_pos == std::wstring::npos || backslash_pos >= source_name.size() - 1)
                    backslash_pos = 0;
                else
                    ++backslash_pos;
                std::size_t dot_pos = source_name.rfind(L'.');
                if (dot_pos == std::wstring::npos || dot_pos < backslash_pos)
                    dot_pos = source_name.size();
                source_name = source_name.substr(backslash_pos, dot_pos - backslash_pos);
            }

            return source_name;
        }

        static LSTATUS create_key(
            HKEY hKey,
            const wchar_t* lpSubKey,
            DWORD Reserved,
            wchar_t* lpClass,
            DWORD dwOptions,
            REGSAM samDesired,
            LPSECURITY_ATTRIBUTES lpSecurityAttributes,
            PHKEY phkResult,
            LPDWORD lpdwDisposition)
        {
            return RegCreateKeyExW(hKey, lpSubKey, Reserved, lpClass, dwOptions, samDesired, lpSecurityAttributes, phkResult, lpdwDisposition);
        }

        static LSTATUS set_value(
            HKEY hKey,
            const wchar_t* lpValueName,
            DWORD Reserved,
            DWORD dwType,
            const BYTE *lpData,
            DWORD cbData)
        {
            return RegSetValueExW(hKey, lpValueName, Reserved, dwType, lpData, cbData);
        }

        static const wchar_t* get_event_message_file_param_name() { return L"EventMessageFile"; }
        static const wchar_t* get_category_message_file_param_name() { return L"CategoryMessageFile"; }
        static const wchar_t* get_category_count_param_name() { return L"CategoryCount"; }
        static const wchar_t* get_types_supported_param_name() { return L"TypesSupported"; }

    };
#endif // BOOST_LOG_USE_WCHAR_T

    //! The structure with parameters that have to be registered in the event log registry key
    template< typename CharT >
    struct registry_params
    {
        typedef std::basic_string< CharT > string_type;

        optional< string_type > event_message_file;
        optional< string_type > category_message_file;
        optional< DWORD > category_count;
        optional< DWORD > types_supported;
    };

    //! A simple guard that closes the registry key on destruction
    struct auto_hkey_close
    {
        explicit auto_hkey_close(HKEY hk) : hk_(hk) {}
        ~auto_hkey_close() { RegCloseKey(hk_); }

    private:
        HKEY hk_;
    };

    //! The function initializes the event log registry key
    template< typename CharT >
    void init_event_log_registry(
        std::basic_string< CharT > const& log_name,
        std::basic_string< CharT > const& source_name,
        bool force,
        registry_params< CharT > const& params)
    {
        typedef std::basic_string< CharT > string_type;
        typedef registry_traits< CharT > registry;
        // Registry key name that contains log description
        string_type reg_key = registry::make_event_log_key(log_name, source_name);

        // Create or open the key
        HKEY hkey = 0;
        DWORD disposition = 0;
        LSTATUS res = registry::create_key(
            HKEY_LOCAL_MACHINE,
            reg_key.c_str(),
            0,
            NULL,
            REG_OPTION_NON_VOLATILE,
            KEY_WRITE,
            NULL,
            &hkey,
            &disposition);
        if (res != ERROR_SUCCESS)
            boost::log::aux::throw_exception(std::runtime_error("Could not create registry key for the event log"));

        auto_hkey_close hkey_guard(hkey);

        if (disposition != REG_OPENED_EXISTING_KEY || force)
        {
            // Fill registry values
            if (!!params.event_message_file)
            {
                // Set the module file name that contains event resources
                string_type const& module_name = params.event_message_file.get();
                res = registry::set_value(
                    hkey,
                    registry::get_event_message_file_param_name(),
                    0,
                    REG_EXPAND_SZ,
                    reinterpret_cast< LPBYTE >(const_cast< CharT* >(module_name.c_str())),
                    static_cast< DWORD >((module_name.size() + 1) * sizeof(CharT)));
                if (res != ERROR_SUCCESS)
                {
                    boost::log::aux::throw_exception(std::runtime_error("Could not create registry value "
                        + log::aux::to_narrow(string_type(registry::get_event_message_file_param_name()))));
                }
            }

            if (!!params.category_message_file)
            {
                // Set the module file name that contains event category resources
                string_type const& module_name = params.category_message_file.get();
                res = registry::set_value(
                    hkey,
                    registry::get_category_message_file_param_name(),
                    0,
                    REG_SZ,
                    reinterpret_cast< LPBYTE >(const_cast< CharT* >(module_name.c_str())),
                    static_cast< DWORD >((module_name.size() + 1) * sizeof(CharT)));
                if (res != ERROR_SUCCESS)
                {
                    boost::log::aux::throw_exception(std::runtime_error("Could not create registry value "
                        + log::aux::to_narrow(string_type(registry::get_category_message_file_param_name()))));
                }
            }

            if (!!params.category_count)
            {
                // Set number of categories
                DWORD category_count = params.category_count.get();
                res = registry::set_value(
                    hkey,
                    registry::get_category_count_param_name(),
                    0,
                    REG_DWORD,
                    reinterpret_cast< LPBYTE >(&category_count),
                    static_cast< DWORD >(sizeof(category_count)));
                if (res != ERROR_SUCCESS)
                {
                    boost::log::aux::throw_exception(std::runtime_error("Could not create registry value "
                        + log::aux::to_narrow(string_type(registry::get_category_count_param_name()))));
                }
            }

            if (!!params.types_supported)
            {
                // Set the supported event types
                DWORD event_types = params.types_supported.get();
                res = registry::set_value(
                    hkey,
                    registry::get_types_supported_param_name(),
                    0,
                    REG_DWORD,
                    reinterpret_cast< LPBYTE >(&event_types),
                    static_cast< DWORD >(sizeof(event_types)));
                if (res != ERROR_SUCCESS)
                {
                    boost::log::aux::throw_exception(std::runtime_error("Could not create registry value "
                        + log::aux::to_narrow(string_type(registry::get_types_supported_param_name()))));
                }
            }
        }
    }

} // namespace aux

} // namespace sinks

} // namespace log

} // namespace boost

#endif // BOOST_LOG_EVENT_LOG_REGISTRY_HPP_INCLUDED_
