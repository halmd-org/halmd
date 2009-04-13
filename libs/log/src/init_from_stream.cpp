/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   init_from_stream.cpp
 * \author Andrey Semashev
 * \date   22.03.2008
 *
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#ifndef BOOST_LOG_NO_SETTINGS_PARSERS_SUPPORT

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <locale>
#include <memory>
#include <utility>
#include <stdexcept>
#include <algorithm>

#if !defined(BOOST_LOG_NO_THREADS) && !defined(BOOST_SPIRIT_THREADSAFE)
#define BOOST_SPIRIT_THREADSAFE
#endif // !defined(BOOST_LOG_NO_THREADS) && !defined(BOOST_SPIRIT_THREADSAFE)

#include <boost/cstdint.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/detail/new_shared.hpp>
#include <boost/log/detail/code_conversion.hpp>
#include <boost/log/detail/throw_exception.hpp>
#include <boost/log/core.hpp>
#include <boost/log/sinks/sink.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sinks/syslog_backend.hpp>
#ifdef BOOST_WINDOWS
#include <boost/log/sinks/debug_output_backend.hpp>
#include <boost/log/sinks/event_log_backend.hpp>
#ifdef BOOST_LOG_USE_WINNT6_API
#include <boost/log/sinks/nt6_event_log_backend.hpp>
#endif // BOOST_LOG_USE_WINNT6_API
#endif // BOOST_WINDOWS
#include <boost/log/detail/singleton.hpp>
#include <boost/log/detail/code_conversion.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/log/utility/rotating_ofstream.hpp>
#include <boost/log/utility/init/from_stream.hpp>
#include <boost/log/utility/init/filter_parser.hpp>
#include <boost/log/utility/init/formatter_parser.hpp>
#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/thread/locks.hpp>
#include <boost/log/detail/light_rw_mutex.hpp>
#include <boost/log/detail/shared_lock_guard.hpp>
#endif
#include "parser_utils.hpp"

#if defined(BOOST_WINDOWS) && defined(BOOST_LOG_USE_WINNT6_API)
#include <guiddef.h>
#endif // defined(BOOST_WINDOWS) && defined(BOOST_LOG_USE_WINNT6_API)

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace {

#if defined(BOOST_WINDOWS) && defined(BOOST_LOG_USE_WINNT6_API)

//! The function parses GUIDS from strings in form "{631EF147-9A84-4e33-8ADD-BCE174C4DBD9}"
template< typename CharT >
GUID parse_guid(std::basic_string< CharT > const& str)
{
    typedef CharT char_type;
    const char_type
        open_bracket = static_cast< char_type >('{'),
        close_bracket = static_cast< char_type >('}'),
        dash = static_cast< char_type >('-');

    const spirit::classic::uint_parser< unsigned char, 16, 2, 2 > octet_p;
    const spirit::classic::uint_parser< unsigned short, 16, 4, 4 > word_p;
    const spirit::classic::uint_parser< unsigned long, 16, 8, 8 > dword_p;

    GUID g;
    spirit::classic::parse_info< const char_type* > result =
        spirit::classic::parse(
            str.c_str(),
            str.c_str() + str.size(),
            (
                spirit::classic::ch_p(open_bracket) >>
                dword_p[spirit::classic::assign_a(g.Data1)] >>
                spirit::classic::ch_p(dash) >>
                word_p[spirit::classic::assign_a(g.Data2)] >>
                spirit::classic::ch_p(dash) >>
                word_p[spirit::classic::assign_a(g.Data3)] >>
                spirit::classic::ch_p(dash) >>
                octet_p[spirit::classic::assign_a(g.Data4[0])] >>
                octet_p[spirit::classic::assign_a(g.Data4[1])] >>
                spirit::classic::ch_p(dash) >>
                octet_p[spirit::classic::assign_a(g.Data4[2])] >>
                octet_p[spirit::classic::assign_a(g.Data4[3])] >>
                octet_p[spirit::classic::assign_a(g.Data4[4])] >>
                octet_p[spirit::classic::assign_a(g.Data4[5])] >>
                octet_p[spirit::classic::assign_a(g.Data4[6])] >>
                octet_p[spirit::classic::assign_a(g.Data4[7])] >>
                spirit::classic::ch_p(close_bracket)
            )
        );

    if (!result.full)
    {
        boost::log::aux::throw_exception(std::runtime_error("Could not recognize CLSID from string "
            + log::aux::to_narrow(str)));
    }

    return g;
}

#endif // defined(BOOST_WINDOWS) && defined(BOOST_LOG_USE_WINNT6_API)

//! The class represents parsed logging settings
template< typename CharT >
class settings
{
public:
    //! Character type
    typedef CharT char_type;
    //! String type
    typedef std::basic_string< char_type > string_type;
    //! The type of the map of parameters and their names
    typedef std::map< string_type, string_type > params_t;
    //! The type of the map of sections
    typedef std::map< string_type, params_t > sections_t;
    //! Structure with character constants
    typedef aux::char_constants< char_type > constants;

private:
    //! Parameters
    sections_t m_Sections;

public:
    //! The constructor reads parameters from the stream
    explicit settings(std::basic_istream< char_type >& strm)
    {
        typedef typename string_type::iterator str_iterator;
        std::locale loc = strm.getloc();
        typename sections_t::iterator current_section = m_Sections.end();
        string_type line;
        for (unsigned int line_counter = 1; strm.good(); ++line_counter)
        {
            line.clear();
            std::getline(strm, line);
            boost::algorithm::trim(line, loc);

            // Skipping empty lines and comments
            // NOTE: The comments are only allowed to be the whole line.
            //       Comments beginning in the middle of the line are not supported.
            if (!line.empty() && line[0] != constants::char_comment)
            {
                // Check if the line is a section starter
                if (line[0] == constants::char_section_bracket_left)
                {
                    str_iterator it = std::find(line.begin() + 1, line.end(), constants::char_section_bracket_right);
                    string_type section_name(line.begin() + 1, it);
                    boost::algorithm::trim(section_name, loc);
                    if (it != line.end() && !section_name.empty())
                    {
                        // Creating a new section
                        current_section = m_Sections.insert(std::make_pair(section_name, params_t())).first;
                    }
                    else
                    {
                        // The section starter is broken
                        std::ostringstream descr;
                        descr << "At line " << line_counter << ". The section header is invalid.";
                        boost::log::aux::throw_exception(std::runtime_error(descr.str()));
                    }
                }
                else
                {
                    // Check that we've already started a section
                    if (current_section != m_Sections.end())
                    {
                        // Find the '=' between the parameter name and value
                        str_iterator it = std::find(line.begin(), line.end(), constants::char_equal);
                        string_type param_name(line.begin(), it);
                        boost::algorithm::trim_right(param_name, loc);
                        if (it != line.end() && !param_name.empty())
                        {
                            // Put the parameter value into the map
                            string_type param_value(++it, line.end());
                            boost::algorithm::trim_left(param_value, loc);
                            if (param_value.size() >= 2
                                && param_value[0] == constants::char_quote && *param_value.rbegin() == constants::char_quote)
                            {
                                param_value = param_value.substr(1, param_value.size() - 2);
                            }

                            current_section->second[param_name] = param_value;
                        }
                        else
                        {
                            // The parameter name is not valid
                            std::ostringstream descr;
                            descr << "At line " << line_counter << ". Parameter description is not valid.";
                            boost::log::aux::throw_exception(std::runtime_error(descr.str()));
                        }
                    }
                    else
                    {
                        // The parameter encountered before any section starter
                        std::ostringstream descr;
                        descr << "At line " << line_counter << ". Parameters are only allowed in sections.";
                        boost::log::aux::throw_exception(std::runtime_error(descr.str()));
                    }
                }
            }
        }
    }

    //! Accessor for the map of sections
    sections_t const& sections() const { return m_Sections; }
};

#ifndef BOOST_FILESYSTEM_NARROW_ONLY

//! A helper trait to generate the appropriate Boost.Filesystem path type
template< typename > struct make_filesystem_path;
template< > struct make_filesystem_path< char > { typedef filesystem::path type; };
template< > struct make_filesystem_path< wchar_t > { typedef filesystem::wpath type; };

#endif // BOOST_FILESYSTEM_NARROW_ONLY

//! The supported sinks repository
template< typename CharT >
struct sinks_repository :
    public log::aux::lazy_singleton< sinks_repository< CharT > >
{
    typedef log::aux::lazy_singleton< sinks_repository< CharT > > base_type;

#if !defined(BOOST_MSVC) || _MSC_VER > 1310
    friend class log::aux::lazy_singleton< sinks_repository< CharT > >;
#else
    friend class base_type;
#endif

    typedef CharT char_type;
    typedef std::basic_string< char_type > string_type;
    typedef boost::log::aux::char_constants< char_type > constants;
    typedef std::map< string_type, string_type > params_t;
    typedef function1< shared_ptr< sinks::sink< char_type > >, params_t const& > sink_factory;
    typedef std::map< string_type, sink_factory > sink_factories;

#if !defined(BOOST_LOG_NO_THREADS)
    //! Synchronization mutex
    log::aux::light_rw_mutex m_Mutex;
#endif
    //! Map of the sink factories
    sink_factories m_Factories;

    //! The function constructs a sink from the settings
    shared_ptr< sinks::sink< char_type > > construct_sink_from_settings(params_t const& params)
    {
        typename params_t::const_iterator dest = params.find(constants::sink_destination_param_name());
        if (dest != params.end())
        {
#if !defined(BOOST_LOG_NO_THREADS)
            log::aux::shared_lock_guard< log::aux::light_rw_mutex > _(m_Mutex);
#endif
            typename sink_factories::const_iterator it = m_Factories.find(dest->second);
            if (it != m_Factories.end())
            {
                return it->second(params);
            }
            else
            {
                boost::log::aux::throw_exception(std::runtime_error("The sink destination is not supported"));
            }
        }
        else
        {
            boost::log::aux::throw_exception(std::runtime_error("The sink destination is not set"));
        }
        // To silence compiler warnings. This return never gets executed.
        return shared_ptr< sinks::sink< char_type > >();
    }

    static void init_instance()
    {
        sinks_repository& instance = base_type::get_instance();
        instance.m_Factories[constants::text_file_destination()] =
            &sinks_repository< char_type >::default_text_file_sink_factory;
        instance.m_Factories[constants::console_destination()] =
            &sinks_repository< char_type >::default_console_sink_factory;
        instance.m_Factories[constants::syslog_destination()] =
            &sinks_repository< char_type >::default_syslog_sink_factory;
#ifdef BOOST_WINDOWS
        instance.m_Factories[constants::debugger_destination()] =
            &sinks_repository< char_type >::default_debugger_sink_factory;
        instance.m_Factories[constants::simple_event_log_destination()] =
            &sinks_repository< char_type >::default_simple_event_log_sink_factory;
#ifdef BOOST_LOG_USE_WINNT6_API
        instance.m_Factories[constants::simple_nt6_event_log_destination()] =
            &sinks_repository< char_type >::default_simple_nt6_event_log_sink_factory;
#endif // BOOST_LOG_USE_WINNT6_API
#endif // BOOST_WINDOWS
    }

private:
    sinks_repository() {}

    //! The function constructs a sink that writes log records to a text file
    static shared_ptr< sinks::sink< char_type > > default_text_file_sink_factory(params_t const& params)
    {
        typedef std::basic_istringstream< char_type > isstream;
        typedef sinks::basic_text_ostream_backend< char_type > backend_t;
        shared_ptr< backend_t > backend = log::aux::new_shared< backend_t >();

        // FileName
#ifndef BOOST_FILESYSTEM_NARROW_ONLY
        typename make_filesystem_path< char_type >::type file_name;
#else
        filesystem::path file_name;
#endif // BOOST_FILESYSTEM_NARROW_ONLY
        typename params_t::const_iterator it = params.find(constants::file_name_param_name());
        if (it != params.end())
        {
#ifndef BOOST_FILESYSTEM_NARROW_ONLY
            file_name = it->second;
#else
            file_name = log::aux::to_narrow(it->second);
#endif // BOOST_FILESYSTEM_NARROW_ONLY
        }
        else
            boost::log::aux::throw_exception(std::runtime_error("File name is not specified"));

        // File rotation params
        shared_ptr< typename backend_t::stream_type > file_stream;
        typename params_t::const_iterator
            rotation_size_param = params.find(constants::rotation_size_param_name()),
            rotation_interval_param = params.find(constants::rotation_interval_param_name());
        unsigned int cond =
            (static_cast< unsigned int >(rotation_size_param != params.end()) << 1)
            | static_cast< unsigned int >(rotation_interval_param != params.end());

        switch (cond)
        {
            case 1:
            {
                // Only rotation interval is set
                unsigned int interval = 0;
                isstream strm(rotation_interval_param->second);
                strm >> interval;

#ifndef BOOST_LOG_BROKEN_STL_ALIGNMENT
                file_stream = log::aux::new_shared< basic_rotating_ofstream< char_type > >(
                    file_name, keywords::rotation_interval = interval);
#else
                file_stream.reset(new basic_rotating_ofstream< char_type > (
                    file_name, keywords::rotation_interval = interval));
#endif // BOOST_LOG_BROKEN_STL_ALIGNMENT

            }
            break;

            case 2:
            {
                // Only rotation size is set
                uintmax_t size = ~static_cast< uintmax_t >(0);
                isstream strm(rotation_size_param->second);
                strm >> size;

#ifndef BOOST_LOG_BROKEN_STL_ALIGNMENT
                file_stream = log::aux::new_shared< basic_rotating_ofstream< char_type > >(
                    file_name, keywords::rotation_size = size);
#else
                file_stream.reset(new basic_rotating_ofstream< char_type >(
                    file_name, keywords::rotation_size = size));
#endif // BOOST_LOG_BROKEN_STL_ALIGNMENT
            }
            break;

            case 3:
            {
                // Both rotation interval and size are set
                unsigned int interval = 0;
                isstream strm_interval(rotation_interval_param->second);
                strm_interval >> interval;

                uintmax_t size = ~static_cast< uintmax_t >(0);
                isstream strm_size(rotation_size_param->second);
                strm_size >> size;

#ifndef BOOST_LOG_BROKEN_STL_ALIGNMENT
                file_stream = log::aux::new_shared< basic_rotating_ofstream< char_type > >(
                    file_name, keywords::rotation_interval = interval, keywords::rotation_size = size);
#else
                file_stream.reset(new basic_rotating_ofstream< char_type >(
                    file_name, keywords::rotation_interval = interval, keywords::rotation_size = size));
#endif // BOOST_LOG_BROKEN_STL_ALIGNMENT
            }
            break;

            default:
            {
                // No rotation required, we can use a simple stream
                typedef filesystem::basic_ofstream< char_type > stream_t;
                shared_ptr< stream_t > p = log::aux::new_shared< stream_t >(
                    file_name, std::ios_base::out | std::ios_base::trunc);
                if (!p->is_open())
                    boost::log::aux::throw_exception(std::runtime_error("Failed to open the destination file"));
                file_stream = p;
            }
        }

        backend->add_stream(file_stream);

        return init_text_ostream_sink(backend, params);
    }

    //! The function constructs a sink that writes log records to the console
    static shared_ptr< sinks::sink< char_type > > default_console_sink_factory(params_t const& params)
    {
        // Construct the backend
        typedef sinks::basic_text_ostream_backend< char_type > backend_t;
        shared_ptr< backend_t > backend = log::aux::new_shared< backend_t >();
        backend->add_stream(
            shared_ptr< typename backend_t::stream_type >(&constants::get_console_log_stream(), empty_deleter()));

        return init_text_ostream_sink(backend, params);
    }

    //! The function constructs a sink that writes log records to the syslog service
    static shared_ptr< sinks::sink< char_type > > default_syslog_sink_factory(params_t const& params)
    {
        // Construct the backend
        typedef sinks::basic_syslog_backend< char_type > backend_t;
        shared_ptr< backend_t > backend = log::aux::new_shared< backend_t >();

        // For now we use only the default level mapping. Will add support for configuration later.
        backend->set_severity_mapper(
            sinks::syslog::basic_direct_severity_mapping< char_type >(constants::default_level_attribute_name()));

        // Setup local and remote addresses
        typename params_t::const_iterator it = params.find(constants::local_address_param_name());
        if (it != params.end())
            backend->set_local_address(log::aux::to_narrow(it->second));

        it = params.find(constants::target_address_param_name());
        if (it != params.end())
            backend->set_target_address(log::aux::to_narrow(it->second));

        return init_sink(backend, params);
    }

#ifdef BOOST_WINDOWS

    //! The function constructs a sink that writes log records to the debugger
    static shared_ptr< sinks::sink< char_type > > default_debugger_sink_factory(params_t const& params)
    {
        // Construct the backend
        typedef sinks::basic_debug_output_backend< char_type > backend_t;
        shared_ptr< backend_t > backend = log::aux::new_shared< backend_t >();

        return init_sink(backend, params);
    }

    //! The function constructs a sink that writes log records to the Windows NT Event Log
    static shared_ptr< sinks::sink< char_type > > default_simple_event_log_sink_factory(params_t const& params)
    {
        typedef sinks::basic_simple_event_log_backend< char_type > backend_t;

        // Determine the log name
        string_type log_name = backend_t::get_default_log_name();
        typename params_t::const_iterator it = params.find(constants::log_name_param_name());
        if (it != params.end())
            log_name = it->second;

        // Determine the log source name
        string_type source_name = backend_t::get_default_source_name();
        it = params.find(constants::source_name_param_name());
        if (it != params.end())
            source_name = it->second;

        // Determine the force flag
        sinks::event_log::registration_mode reg_mode = sinks::event_log::on_demand;
        it = params.find(constants::registration_param_name());
        if (it != params.end())
        {
            if (it->second == constants::registration_never())
                reg_mode = sinks::event_log::never;
            else if (it->second == constants::registration_on_demand())
                reg_mode = sinks::event_log::on_demand;
            else if (it->second == constants::registration_forced())
                reg_mode = sinks::event_log::forced;
            else
            {
                boost::log::aux::throw_exception(std::runtime_error(
                    "The registration mode \"" + log::aux::to_narrow(it->second) + "\" is not supported"));
            }
        }

        // Construct the backend
        shared_ptr< backend_t > backend(new backend_t((
            keywords::log_name = log_name,
            keywords::log_source = source_name,
            keywords::registration = reg_mode)));

        // For now we use only the default event type mapping. Will add support for configuration later.
        backend->set_event_type_mapper(
            sinks::event_log::basic_direct_event_type_mapping< char_type >(constants::default_level_attribute_name()));

        return init_sink(backend, params);
    }

#ifdef BOOST_LOG_USE_WINNT6_API

    //! The function constructs a sink that writes log records to the Windows NT 6 Event Log
    static shared_ptr< sinks::sink< char_type > > default_simple_nt6_event_log_sink_factory(params_t const& params)
    {
        typedef experimental::sinks::basic_simple_nt6_event_log_backend< char_type > backend_t;

        // Determine the provider GUID
        GUID provider_id = backend_t::get_default_provider_id();
        typename params_t::const_iterator it = params.find(constants::provider_id_param_name());
        if (it != params.end())
            provider_id = parse_guid(it->second);

        // Construct the backend
        shared_ptr< backend_t > backend(log::aux::new_shared< backend_t >(boost::cref(provider_id)));

        // For now we use only the default level mapping. Will add support for configuration later.
        backend->set_severity_mapper(
            experimental::sinks::etw::basic_direct_severity_mapping< char_type >(constants::default_level_attribute_name()));

        return init_sink(backend, params);
    }

#endif // BOOST_LOG_USE_WINNT6_API
#endif // BOOST_WINDOWS

    //! The function initializes common parameters of text stream sink and returns the constructed sink
    static shared_ptr< sinks::sink< char_type > > init_text_ostream_sink(
        shared_ptr< sinks::basic_text_ostream_backend< char_type > > const& backend, params_t const& params)
    {
        typedef sinks::basic_text_ostream_backend< char_type > backend_t;

        // AutoFlush
        typedef std::basic_istringstream< char_type > isstream;
        typename params_t::const_iterator it = params.find(constants::auto_flush_param_name());
        if (it != params.end())
        {
            isstream strm(it->second);
            strm.setf(std::ios_base::boolalpha);
            bool f = false;
            strm >> f;
            backend->auto_flush(f);
        }

        return init_sink(backend, params);
    }

    //! The function initializes common parameters of a formatting sink and returns the constructed sink
    template< typename BackendT >
    static shared_ptr< sinks::sink< char_type > > init_sink(shared_ptr< BackendT > const& backend, params_t const& params)
    {
        typedef BackendT backend_t;
        typedef std::basic_istringstream< char_type > isstream;

        // Filter
        typedef typename sinks::sink< char_type >::filter_type filter_type;
        filter_type filt;
        typename params_t::const_iterator it = params.find(constants::filter_param_name());
        if (it != params.end())
        {
            filt = parse_filter(it->second);
        }

        // Formatter
        it = params.find(constants::format_param_name());
        if (it != params.end())
        {
            backend->set_formatter(parse_formatter(it->second));
        }

        shared_ptr< sinks::sink< char_type > > p;

#if !defined(BOOST_LOG_NO_THREADS)
        // Asynchronous
        bool async = false;
        it = params.find(constants::asynchronous_param_name());
        if (it != params.end())
        {
            isstream strm(it->second);
            strm.setf(std::ios_base::boolalpha);
            strm >> async;
        }

        // Construct the frontend, considering Asynchronous parameter
        if (!async)
            p = log::aux::new_shared< sinks::synchronous_sink< backend_t > >(backend);
        else
            p = log::aux::new_shared< sinks::asynchronous_sink< backend_t > >(backend);
#else
        // When multithreading is disabled we always use the unlocked sink frontend
        p = log::aux::new_shared< sinks::unlocked_sink< backend_t > >(backend);
#endif

        p->set_filter(filt);

        return p;
    }
};

//! The function applies the settings to the logging core
template< typename CharT >
void apply_core_settings(std::map< std::basic_string< CharT >, std::basic_string< CharT > > const& params)
{
    typedef CharT char_type;
    typedef std::basic_string< char_type > string_type;
    typedef std::map< string_type, string_type > params_t;
    typedef aux::char_constants< char_type > constants;
    typedef std::basic_istringstream< char_type > isstream;
    typedef basic_core< char_type > core_t;
    shared_ptr< core_t > core = core_t::get();

    // Filter
    typename params_t::const_iterator it =
        params.find(constants::filter_param_name());
    if (it != params.end())
        core->set_filter(parse_filter(it->second));
    else
        core->reset_filter();

    // DisableLogging
    it = params.find(constants::core_disable_logging_param_name());
    if (it != params.end())
    {
        isstream strm(it->second);
        strm.setf(std::ios_base::boolalpha);
        bool f = false;
        strm >> f;
        core->set_logging_enabled(!f);
    }
    else
        core->set_logging_enabled(true);
}

} // namespace

//! The function registers a factory for a sink
template< typename CharT >
void register_sink_factory(
    const CharT* sink_name,
    function1<
        shared_ptr< sinks::sink< CharT > >,
        std::map< std::basic_string< CharT >, std::basic_string< CharT > > const&
    > const& factory)
{
    sinks_repository< CharT >& repo = sinks_repository< CharT >::get();
#if !defined(BOOST_LOG_NO_THREADS)
    lock_guard< log::aux::light_rw_mutex > _(repo.m_Mutex);
#endif
    repo.m_Factories[sink_name] = factory;
}


//! The function initializes the logging library from a stream containing logging settings
template< typename CharT >
void init_from_stream(std::basic_istream< CharT >& strm)
{
    typedef CharT char_type;
    typedef std::basic_string< char_type > string_type;
    typedef basic_core< char_type > core_t;
    typedef sinks_repository< char_type > sinks_repo_t;

    // Parse the settings
    typedef settings< CharT > settings_t;
    typedef typename settings_t::constants constants;
    settings_t setts(strm);

    // Apply core settings
    typename settings_t::sections_t const& sections = setts.sections();
    typename settings_t::sections_t::const_iterator it =
        sections.find(constants::core_section_name());
    if (it != sections.end())
        apply_core_settings(it->second);

    // Construct and initialize sinks
    sinks_repo_t& sinks_repo = sinks_repo_t::get();
    string_type sink_prefix = constants::sink_section_name_prefix();
    std::vector< shared_ptr< sinks::sink< char_type > > > new_sinks;
    for (it = setts.sections().begin(); it != setts.sections().end(); ++it)
    {
        if (it->first.compare(0, sink_prefix.size(), sink_prefix) == 0)
            new_sinks.push_back(sinks_repo.construct_sink_from_settings(it->second));
    }
    std::for_each(new_sinks.begin(), new_sinks.end(), bind(&core_t::add_sink, core_t::get(), _1));
}


#ifdef BOOST_LOG_USE_CHAR
template BOOST_LOG_EXPORT
void register_sink_factory< char >(
    const char* sink_name,
    function1<
        shared_ptr< sinks::sink< char > >,
        std::map< std::basic_string< char >, std::basic_string< char > > const&
    > const& factory);
template BOOST_LOG_EXPORT void init_from_stream< char >(std::basic_istream< char >& strm);
#endif

#ifdef BOOST_LOG_USE_WCHAR_T
template BOOST_LOG_EXPORT
void register_sink_factory< wchar_t >(
    const wchar_t* sink_name,
    function1<
        shared_ptr< sinks::sink< wchar_t > >,
        std::map< std::basic_string< wchar_t >, std::basic_string< wchar_t > > const&
    > const& factory);
template BOOST_LOG_EXPORT void init_from_stream< wchar_t >(std::basic_istream< wchar_t >& strm);
#endif

} // namespace log

} // namespace boost

#endif // BOOST_LOG_NO_SETTINGS_PARSERS_SUPPORT
