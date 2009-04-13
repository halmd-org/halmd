/*!
 * (C) 2008 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   event_log_backend.cpp
 * \author Andrey Semashev
 * \date   07.11.2008
 *
 * \brief  A logging sink backend that uses Windows NT event log API
 *         for signalling application events.
 */

#define WIN32_LEAN_AND_MEAN

#include "windows_version.hpp"
#include <windows.h>
#include <psapi.h>
#include "simple_event_log.h"

#ifdef _MSC_VER
#pragma comment(lib, "psapi.lib")
#pragma comment(lib, "advapi32.lib")
#endif

#include <memory>
#include <string>
#include <vector>
#include <ostream>
#include <stdexcept>
#include <boost/scoped_array.hpp>
#include <boost/log/sinks/event_log_backend.hpp>
#include <boost/log/sinks/event_log_constants.hpp>
#include <boost/log/detail/throw_exception.hpp>
#include <boost/log/detail/cleanup_scope_guard.hpp>
#include <boost/log/detail/attachable_sstream_buf.hpp>
#ifndef BOOST_LOG_NO_THREADS
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/thread/once.hpp>
#endif // BOOST_LOG_NO_THREADS
#include "event_log_registry.hpp"

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sinks {

namespace event_log {

    //  Windows event types
    BOOST_LOG_EXPORT const event_type_t success = { EVENTLOG_SUCCESS };
    BOOST_LOG_EXPORT const event_type_t info = { EVENTLOG_INFORMATION_TYPE };
    BOOST_LOG_EXPORT const event_type_t warning = { EVENTLOG_WARNING_TYPE };
    BOOST_LOG_EXPORT const event_type_t error = { EVENTLOG_ERROR_TYPE };

} // namespace event_log

namespace {

#ifdef BOOST_LOG_USE_CHAR
    //! A simple forwarder to the ReportEvent API
    inline BOOL report_event(
        HANDLE hEventLog,
        WORD wType,
        WORD wCategory,
        DWORD dwEventID,
        PSID lpUserSid,
        WORD wNumStrings,
        DWORD dwDataSize,
        const char** lpStrings,
        LPVOID lpRawData)
    {
        return ReportEventA(hEventLog, wType, wCategory, dwEventID, lpUserSid, wNumStrings, dwDataSize, lpStrings, lpRawData);
    }
    //! A simple forwarder to the GetModuleFileName API
    inline DWORD get_module_file_name(HMODULE hModule, char* lpFilename, DWORD nSize)
    {
        return GetModuleFileNameA(hModule, lpFilename, nSize);
    }
    //! A simple forwarder to the RegisterEventSource API
    inline HANDLE register_event_source(const char* lpUNCServerName, const char* lpSourceName)
    {
        return RegisterEventSourceA(lpUNCServerName, lpSourceName);
    }
    //! The function completes default source name for the sink backend
    inline void complete_default_simple_event_log_source_name(std::string& name)
    {
        name += " simple event source";
    }
    //! The function completes default source name for the sink backend
    inline void complete_default_event_log_source_name(std::string& name)
    {
        name += " event source";
    }
#endif // BOOST_LOG_USE_CHAR

#ifdef BOOST_LOG_USE_WCHAR_T
    //! A simple forwarder to the ReportEvent API
    inline BOOL report_event(
        HANDLE hEventLog,
        WORD wType,
        WORD wCategory,
        DWORD dwEventID,
        PSID lpUserSid,
        WORD wNumStrings,
        DWORD dwDataSize,
        const wchar_t** lpStrings,
        LPVOID lpRawData)
    {
        return ReportEventW(hEventLog, wType, wCategory, dwEventID, lpUserSid, wNumStrings, dwDataSize, lpStrings, lpRawData);
    }
    //! A simple forwarder to the GetModuleFileName API
    inline DWORD get_module_file_name(HMODULE hModule, wchar_t* lpFilename, DWORD nSize)
    {
        return GetModuleFileNameW(hModule, lpFilename, nSize);
    }
    //! A simple forwarder to the RegisterEventSource API
    inline HANDLE register_event_source(const wchar_t* lpUNCServerName, const wchar_t* lpSourceName)
    {
        return RegisterEventSourceW(lpUNCServerName, lpSourceName);
    }
    //! The function completes default source name for the sink backend
    inline void complete_default_simple_event_log_source_name(std::wstring& name)
    {
        name += L" simple event source";
    }
    //! The function completes default source name for the sink backend
    inline void complete_default_event_log_source_name(std::wstring& name)
    {
        name += L" event source";
    }
#endif // BOOST_LOG_USE_WCHAR_T

    //! The function finds the handle for the current module
    void init_self_module_handle(HMODULE& handle)
    {
        // Acquire all modules of the current process
        HANDLE hProcess = GetCurrentProcess();
        std::vector< HMODULE > modules;
        DWORD module_count = 1024;
        do
        {
            modules.resize(module_count, HMODULE(0));
            BOOL res = EnumProcessModules(
                hProcess,
                &modules[0],
                static_cast< DWORD >(modules.size() * sizeof(HMODULE)),
                &module_count);
            module_count /= sizeof(HMODULE);

            if (!res)
                boost::log::aux::throw_exception(std::runtime_error("Could not enumerate process modules"));
        }
        while (module_count > modules.size());
        modules.resize(module_count, HMODULE(0));

        // Now find the current module among them
        void* p = (void*)&init_self_module_handle;
        for (std::size_t i = 0, n = modules.size(); i < n; ++i)
        {
            MODULEINFO info;
            if (!GetModuleInformation(hProcess, modules[i], &info, sizeof(info)))
                boost::log::aux::throw_exception(std::runtime_error("Could not acquire module information"));

            if (info.lpBaseOfDll <= p && (static_cast< unsigned char* >(info.lpBaseOfDll) + info.SizeOfImage) > p)
            {
                // Found it
                handle = modules[i];
                break;
            }
        }

        if (!handle)
            boost::log::aux::throw_exception(std::runtime_error("Could not find self module information"));
    }

    //! Retrieves the full name of the current module (be that dll or exe)
    template< typename CharT >
    std::basic_string< CharT > get_current_module_name()
    {
        static HMODULE hSelfModule = 0;
#ifndef BOOST_LOG_NO_THREADS
        static once_flag flag = BOOST_ONCE_INIT;
        boost::call_once(flag, boost::bind(&init_self_module_handle, boost::ref(hSelfModule)));
#else
        if (!hSelfModule)
            init_self_module_handle(hSelfModule);
#endif // BOOST_LOG_NO_THREADS

        // Get the module file name
        CharT buf[MAX_PATH];
        DWORD size = get_module_file_name(hSelfModule, buf, sizeof(buf) / sizeof(*buf));
        if (size == 0)
            boost::log::aux::throw_exception(std::runtime_error("Could not get module file name"));

        return std::basic_string< CharT >(buf, buf + size);
    }

} // namespace

//////////////////////////////////////////////////////////////////////////
//  Simple event log backend implementation
//////////////////////////////////////////////////////////////////////////
//! Sink backend implementation
template< typename CharT >
struct basic_simple_event_log_backend< CharT >::implementation
{
    //! A handle for the registered event provider
    HANDLE m_SourceHandle;
    //! A level mapping functor
    event_type_mapper_type m_LevelMapper;

    implementation() : m_SourceHandle(0)
    {
    }
};

//! Default constructor. Registers event source Boost.Log <Boost version> in the Application log.
template< typename CharT >
basic_simple_event_log_backend< CharT >::basic_simple_event_log_backend() :
    m_pImpl(construct(string_type(), get_default_log_name(), get_default_source_name(), event_log::on_demand))
{
}

//! Destructor
template< typename CharT >
basic_simple_event_log_backend< CharT >::~basic_simple_event_log_backend()
{
    DeregisterEventSource(m_pImpl->m_SourceHandle);
    delete m_pImpl;
}

//! Constructs backend implementation
template< typename CharT >
typename basic_simple_event_log_backend< CharT >::implementation*
basic_simple_event_log_backend< CharT >::construct(
    string_type const& target, string_type const& log_name, string_type const& source_name, event_log::registration_mode reg_mode)
{
    if (reg_mode != event_log::never)
    {
        aux::registry_params< char_type > reg_params;
        reg_params.event_message_file = get_current_module_name< char_type >();
        reg_params.types_supported = DWORD(
            EVENTLOG_SUCCESS |
            EVENTLOG_INFORMATION_TYPE |
            EVENTLOG_WARNING_TYPE |
            EVENTLOG_ERROR_TYPE);
        aux::init_event_log_registry(log_name, source_name, reg_mode == event_log::forced, reg_params);
    }

    std::auto_ptr< implementation > p(new implementation());

    const char_type* target_unc = NULL;
    if (!target.empty())
        target_unc = target.c_str();

    HANDLE hSource = register_event_source(target_unc, source_name.c_str());
    if (!hSource)
        boost::log::aux::throw_exception(std::runtime_error("Could not register event source"));

    p->m_SourceHandle = hSource;

    return p.release();
}

//! Returns default log name
template< typename CharT >
typename basic_simple_event_log_backend< CharT >::string_type
basic_simple_event_log_backend< CharT >::get_default_log_name()
{
    return aux::registry_traits< char_type >::make_default_log_name();
}

//! Returns default source name
template< typename CharT >
typename basic_simple_event_log_backend< CharT >::string_type
basic_simple_event_log_backend< CharT >::get_default_source_name()
{
    string_type source_name = aux::registry_traits< char_type >::make_default_source_name();
    complete_default_simple_event_log_source_name(source_name);
    return source_name;
}

//! The method installs the function object that maps application severity levels to WinAPI event types
template< typename CharT >
void basic_simple_event_log_backend< CharT >::set_event_type_mapper(event_type_mapper_type const& mapper)
{
    m_pImpl->m_LevelMapper = mapper;
}

//! The method puts the formatted message to the event log
template< typename CharT >
void basic_simple_event_log_backend< CharT >::do_consume(record_type const& record, target_string_type const& formatted_message)
{
    const char_type* message = formatted_message.c_str();
    WORD event_type = EVENTLOG_INFORMATION_TYPE;
    if (!m_pImpl->m_LevelMapper.empty())
        event_type = m_pImpl->m_LevelMapper(record.attribute_values()).value;

    DWORD event_id;
    switch (event_type)
    {
    case EVENTLOG_SUCCESS:
        event_id = BOOST_LOG_MSG_DEBUG; break;
    case EVENTLOG_WARNING_TYPE:
        event_id = BOOST_LOG_MSG_WARNING; break;
    case EVENTLOG_ERROR_TYPE:
        event_id = BOOST_LOG_MSG_ERROR; break;
    case EVENTLOG_INFORMATION_TYPE:
    default:
        event_id = BOOST_LOG_MSG_INFO; break;
    }

    report_event(
        m_pImpl->m_SourceHandle,    // Event log handle.
        event_type,                 // Event type.
        0,                          // Event category.
        event_id,                   // Event identifier.
        NULL,                       // No user security identifier.
        1,                          // Number of substitution strings.
        0,                          // No data.
        &message,                   // Pointer to strings.
        NULL);                      // No data.
}

//////////////////////////////////////////////////////////////////////////
//  Customizable event log backend implementation
//////////////////////////////////////////////////////////////////////////
namespace event_log {

    template< typename CharT >
    class basic_event_composer< CharT >::insertion_composer
    {
    public:
        //! Function object result type
        typedef void result_type;

    private:
        //! The list of insertion composers (in backward order)
        typedef std::vector< formatter_type > formatters;

    private:
        //! The insertion string composers
        formatters m_Formatters;

    public:
        //! Default constructor
        insertion_composer() {}
        //! Composition operator
        void operator() (
            values_view_type const& attributes,
            string_type const& message,
            insertion_list& insertions) const
        {
            std::size_t size = m_Formatters.size();
            insertions.resize(size);
            for (std::size_t i = 0; i < size; ++i)
            {
                log::aux::basic_ostringstreambuf< char_type > buf(insertions[i]);
                stream_type strm(&buf);
                m_Formatters[i](strm, attributes, message);
                strm.flush();
            }
        }
        //! Adds a new formatter to the list
        void add_formatter(formatter_type const& fmt)
        {
            m_Formatters.push_back(formatter_type(fmt));
        }
    };

    //! Default constructor
    template< typename CharT >
    basic_event_composer< CharT >::basic_event_composer(event_id_mapper_type const& id_mapper) :
        m_EventIDMapper(id_mapper)
    {
    }
    //! Copy constructor
    template< typename CharT >
    basic_event_composer< CharT >::basic_event_composer(basic_event_composer const& that) :
        m_EventIDMapper(that.m_EventIDMapper),
        m_EventMap(that.m_EventMap)
    {
    }
    //! Destructor
    template< typename CharT >
    basic_event_composer< CharT >::~basic_event_composer()
    {
    }

    //! Assignment
    template< typename CharT >
    basic_event_composer< CharT >& basic_event_composer< CharT >::operator= (basic_event_composer that)
    {
        swap(that);
        return *this;
    }
    //! Swapping
    template< typename CharT >
    void basic_event_composer< CharT >::swap(basic_event_composer& that)
    {
        m_EventIDMapper.swap(that.m_EventIDMapper);
        m_EventMap.swap(that.m_EventMap);
    }
    //! Creates a new entry for a message
    template< typename CharT >
    typename basic_event_composer< CharT >::event_map_reference
    basic_event_composer< CharT >::operator[] (event_id_t id)
    {
        return event_map_reference(id, *this);
    }
    //! Creates a new entry for a message
    template< typename CharT >
    typename basic_event_composer< CharT >::event_map_reference
    basic_event_composer< CharT >::operator[] (event_id_t::integer_type id)
    {
        return event_map_reference(make_event_id(id), *this);
    }

    //! Event composition operator
    template< typename CharT >
    event_id_t basic_event_composer< CharT >::operator() (
        values_view_type const& attributes,
        string_type const& message,
        insertion_list& insertions) const
    {
        event_id_t id = m_EventIDMapper(attributes);
        typename event_map::const_iterator it = m_EventMap.find(id);
        if (it != m_EventMap.end())
            it->second(attributes, message, insertions);
        return id;
    }

    //! Adds a formatter to the insertion composers list
    template< typename CharT >
    typename basic_event_composer< CharT >::insertion_composer*
    basic_event_composer< CharT >::add_formatter(
        event_id_t id, insertion_composer* composer, formatter_type const& fmt)
    {
        if (!composer)
            composer = &m_EventMap[id];
        composer->add_formatter(fmt);
        return composer;
    }

#ifdef BOOST_LOG_USE_CHAR
    template class BOOST_LOG_EXPORT basic_event_composer< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
    template class BOOST_LOG_EXPORT basic_event_composer< wchar_t >;
#endif

} // namespace event_log


//! Backend implementation
template< typename CharT >
struct basic_event_log_backend< CharT >::implementation
{
    //  NOTE: This order of data members is critical for MSVC 9.0 in debug mode,
    //        as it ICEs if boost::functions are not the first members. Doh!

    //! An event category mapper
    event_category_mapper_type m_CategoryMapper;
    //! A level mapping functor
    event_type_mapper_type m_LevelMapper;

    //! A handle for the registered event provider
    HANDLE m_SourceHandle;
    //! A functor that composes an event
    event_composer_type m_EventComposer;
    //! An array of formatted insertions
    insertion_list m_Insertions;

    implementation() : m_SourceHandle(0)
    {
    }
};

//! Destructor
template< typename CharT >
basic_event_log_backend< CharT >::~basic_event_log_backend()
{
    DeregisterEventSource(m_pImpl->m_SourceHandle);
    delete m_pImpl;
}

//! Constructs backend implementation
template< typename CharT >
typename basic_event_log_backend< CharT >::implementation* basic_event_log_backend< CharT >::construct(
    filesystem::path const& message_file_name,
    string_type const& target,
    string_type const& log_name,
    string_type const& source_name,
    event_log::registration_mode reg_mode)
{
    string_type fname;
    log::aux::code_convert(message_file_name.string(), fname);
    return construct_impl(fname, target, log_name, source_name, reg_mode);
}
#ifndef BOOST_FILESYSTEM_NARROW_ONLY
//! Constructs backend implementation
template< typename CharT >
typename basic_event_log_backend< CharT >::implementation* basic_event_log_backend< CharT >::construct(
    filesystem::wpath const& message_file_name,
    string_type const& target,
    string_type const& log_name,
    string_type const& source_name,
    event_log::registration_mode reg_mode)
{
    string_type fname;
    log::aux::code_convert(message_file_name.string(), fname);
    return construct_impl(fname, target, log_name, source_name, reg_mode);
}
#endif // BOOST_FILESYSTEM_NARROW_ONLY
//! Constructs backend implementation
template< typename CharT >
typename basic_event_log_backend< CharT >::implementation* basic_event_log_backend< CharT >::construct(
    std::string const& message_file_name,
    string_type const& target,
    string_type const& log_name,
    string_type const& source_name,
    event_log::registration_mode reg_mode)
{
    string_type fname;
    log::aux::code_convert(message_file_name, fname);
    return construct_impl(fname, target, log_name, source_name, reg_mode);
}
//! Constructs backend implementation
template< typename CharT >
typename basic_event_log_backend< CharT >::implementation* basic_event_log_backend< CharT >::construct(
    std::wstring const& message_file_name,
    string_type const& target,
    string_type const& log_name,
    string_type const& source_name,
    event_log::registration_mode reg_mode)
{
    string_type fname;
    log::aux::code_convert(message_file_name, fname);
    return construct_impl(fname, target, log_name, source_name, reg_mode);
}

//! Constructs backend implementation
template< typename CharT >
typename basic_event_log_backend< CharT >::implementation* basic_event_log_backend< CharT >::construct_impl(
    string_type const& message_file_name,
    string_type const& target,
    string_type const& log_name,
    string_type const& source_name,
    event_log::registration_mode reg_mode)
{
    if (reg_mode != event_log::never)
    {
        aux::registry_params< char_type > reg_params;
        reg_params.event_message_file = message_file_name;
        reg_params.types_supported = DWORD(
            EVENTLOG_SUCCESS |
            EVENTLOG_INFORMATION_TYPE |
            EVENTLOG_WARNING_TYPE |
            EVENTLOG_ERROR_TYPE);
        aux::init_event_log_registry(log_name, source_name, reg_mode == event_log::forced, reg_params);
    }

    std::auto_ptr< implementation > p(new implementation());

    const char_type* target_unc = NULL;
    if (!target.empty())
        target_unc = target.c_str();

    HANDLE hSource = register_event_source(target_unc, source_name.c_str());
    if (!hSource)
        boost::log::aux::throw_exception(std::runtime_error("Could not register event source"));

    p->m_SourceHandle = hSource;

    return p.release();
}

//! The method puts the formatted message to the event log
template< typename CharT >
void basic_event_log_backend< CharT >::consume(record_type const& record)
{
    if (!m_pImpl->m_EventComposer.empty())
    {
        log::aux::cleanup_guard< insertion_list > _(m_pImpl->m_Insertions);
        values_view_type const& values = record.attribute_values();

        // Get event ID and construct insertions
        DWORD event_id = m_pImpl->m_EventComposer(values, record.message(), m_pImpl->m_Insertions).value;
        WORD string_count = static_cast< WORD >(m_pImpl->m_Insertions.size());
        scoped_array< const char_type* > strings(new const char_type*[string_count]);
        for (WORD i = 0; i < string_count; ++i)
            strings[i] = m_pImpl->m_Insertions[i].c_str();

        // Get event type
        WORD event_type = EVENTLOG_INFORMATION_TYPE;
        if (!m_pImpl->m_LevelMapper.empty())
            event_type = m_pImpl->m_LevelMapper(values).value;

        WORD event_category = 0;
        if (!m_pImpl->m_CategoryMapper.empty())
            event_category = m_pImpl->m_CategoryMapper(values).value;

        report_event(
            m_pImpl->m_SourceHandle,       // Event log handle.
            event_type,                    // Event type.
            event_category,                // Event category.
            event_id,                      // Event identifier.
            NULL,                          // No user security identifier.
            string_count,                  // Number of substitution strings.
            0,                             // No data.
            strings.get(),                 // Pointer to strings.
            NULL);                         // No data.
    }
}

//! Returns default log name
template< typename CharT >
typename basic_event_log_backend< CharT >::string_type
basic_event_log_backend< CharT >::get_default_log_name()
{
    return aux::registry_traits< char_type >::make_default_log_name();
}

//! Returns default source name
template< typename CharT >
typename basic_event_log_backend< CharT >::string_type
basic_event_log_backend< CharT >::get_default_source_name()
{
    string_type source_name = aux::registry_traits< char_type >::make_default_source_name();
    complete_default_event_log_source_name(source_name);
    return source_name;
}

//! The method installs the function object that maps application severity levels to WinAPI event types
template< typename CharT >
void basic_event_log_backend< CharT >::set_event_type_mapper(event_type_mapper_type const& mapper)
{
    m_pImpl->m_LevelMapper = mapper;
}

//! The method installs the function object that extracts event category from attribute values
template< typename CharT >
void basic_event_log_backend< CharT >::set_event_category_mapper(event_category_mapper_type const& mapper)
{
    m_pImpl->m_CategoryMapper = mapper;
}

/*!
 * The method installs the function object that extracts event identifier from the attributes and creates
 * insertion strings that will replace placeholders in the event message.
 */
template< typename CharT >
void basic_event_log_backend< CharT >::set_event_composer(event_composer_type const& composer)
{
    m_pImpl->m_EventComposer = composer;
}


#ifdef BOOST_LOG_USE_CHAR
template class BOOST_LOG_EXPORT basic_simple_event_log_backend< char >;
template class BOOST_LOG_EXPORT basic_event_log_backend< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class BOOST_LOG_EXPORT basic_simple_event_log_backend< wchar_t >;
template class BOOST_LOG_EXPORT basic_event_log_backend< wchar_t >;
#endif

} // namespace sinks

} // namespace log

} // namespace boost
