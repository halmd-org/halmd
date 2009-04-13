/*
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * This header is the Boost.Log library implementation, see the library documentation
 * at http://www.boost.org/libs/log/doc/log.html.
 */
/*!
 * \file   rotating_ofstream.hpp
 * \author Andrey Semashev
 * \date   29.07.2007
 *
 * The header contains implementation of a rotating file stream.
 */

#if (defined(_MSC_VER) && _MSC_VER > 1000)
#pragma once
#endif // _MSC_VER > 1000

#ifndef BOOST_LOG_UTILITY_ROTATING_OFSTREAM_HPP_INCLUDED_
#define BOOST_LOG_UTILITY_ROTATING_OFSTREAM_HPP_INCLUDED_

#include <ios>
#include <memory>
#include <string>
#include <ostream>
#include <stdexcept>
#include <boost/ref.hpp>
#include <boost/limits.hpp>
#include <boost/optional.hpp>
#include <boost/cstdint.hpp>
#include <boost/compatibility/cpp_c_headers/cctype>
#include <boost/compatibility/cpp_c_headers/cwctype>
#include <boost/compatibility/cpp_c_headers/ctime>
#include <boost/compatibility/cpp_c_headers/cstdio>
#include <boost/compatibility/cpp_c_headers/cstdlib>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>
#ifndef BOOST_FILESYSTEM_NARROW_ONLY
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#endif // BOOST_FILESYSTEM_NARROW_ONLY
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/utility/base_from_member.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/function/function1.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/log/detail/prologue.hpp>
#include <boost/log/detail/throw_exception.hpp>
#include <boost/log/detail/cleanup_scope_guard.hpp>
#include <boost/log/detail/code_conversion.hpp>
#include <boost/log/detail/snprintf.hpp>
#include <boost/log/attributes/time_traits.hpp>
#include <boost/log/utility/record_writer.hpp>
#include <boost/log/keywords/rotation_size.hpp>
#include <boost/log/keywords/rotation_interval.hpp>
#include <boost/log/keywords/open_mode.hpp>

#ifdef _MSC_VER
#pragma warning(push)
// 'm_A' : class 'A' needs to have dll-interface to be used by clients of class 'B'
#pragma warning(disable: 4251)
// non dll-interface struct 'A' used as base for dll-interface class 'B'
#pragma warning(disable: 4275)
#endif // _MSC_VER

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace aux {

    //! Stream buffer type generator
    template<
        typename CharT,
        typename TraitsT = std::char_traits< CharT >,
        typename AllocatorT = std::allocator< CharT >
    >
    struct make_narrowing_ostringstreambuf;

    template< typename TraitsT, typename AllocatorT >
    struct make_narrowing_ostringstreambuf< char, TraitsT, AllocatorT >
    {
        typedef basic_ostringstreambuf< char, TraitsT, AllocatorT > type;
    };

    template< typename TraitsT, typename AllocatorT >
    struct make_narrowing_ostringstreambuf< wchar_t, TraitsT, AllocatorT >
    {
        typedef converting_ostringstreambuf< wchar_t, TraitsT > type;
    };


    //! An auxiliary traits that contain various constants and functions regarding string and character operations
    template< typename CharT >
    struct file_controller_traits;

#ifdef BOOST_LOG_USE_CHAR
    template< >
    struct file_controller_traits< char >
    {
        enum
        {
            percent = '%',
            number_placeholder = 'N',
            digit_placeholder = 'u',
            space = ' ',
            plus = '+',
            minus = '-',
            zero = '0',
            dot = '.'
        };
        static bool is_digit(char c)
        {
            return (isdigit(c) != 0);
        }

        typedef int (*lazy_sprintf)(char*, std::size_t, const char*, ...);

        static lazy_sprintf get_lazy_sprintf()
        {
            return (lazy_sprintf)&boost::log::aux::snprintf;
        }
    };
#endif // BOOST_LOG_USE_CHAR

#ifdef BOOST_LOG_USE_WCHAR_T
    template< >
    struct file_controller_traits< wchar_t >
    {
        enum
        {
            percent = L'%',
            number_placeholder = L'N',
            digit_placeholder = L'u',
            space = L' ',
            plus = L'+',
            minus = L'-',
            zero = L'0',
            dot = L'.'
        };
        static bool is_digit(wchar_t c)
        {
            return (iswdigit(c) != 0);
        }

        typedef int (*lazy_sprintf)(wchar_t*, std::size_t, const wchar_t*, ...);

        static lazy_sprintf get_lazy_sprintf()
        {
            return (lazy_sprintf)&boost::log::aux::swprintf;
        }
    };
#endif // BOOST_LOG_USE_WCHAR_T

    //! Base class for file controller instantiations
    struct BOOST_LOG_NO_VTABLE file_controller_base
    {
    public:
        //! The type of the handler of the event of opening of a new file
        typedef function1< void, std::ostream& > open_handler_type;
        //! The type of the handler of the event of closing of a file
        typedef function1< void, std::ostream& > close_handler_type;

    protected:
        //! File buffer
        filesystem::basic_filebuf< char > m_File;
        //! File open mode
        std::ios_base::openmode m_OpenMode;

        //! Chars written to the file
        uintmax_t m_Written;
        //! File size rotation limit
        uintmax_t m_RotationSize;

        //! Last rotation time
        std::time_t m_LastRotation;
        //! Rotation time interval in seconds
        unsigned int m_RotationInterval;

    public:
        file_controller_base(std::ios_base::openmode mode, uintmax_t rot_size, unsigned int rot_int) :
            m_OpenMode(mode),
            m_Written(0),
            m_RotationSize(rot_size),
            m_LastRotation(0),
            m_RotationInterval(rot_int)
        {
            m_OpenMode |= std::ios_base::out | std::ios_base::trunc;
            m_OpenMode &= ~std::ios_base::in;
        }
        virtual ~file_controller_base() {}

        //! The function flushes the buffered data to the file
        int flush()
        {
            return m_File.pubsync();
        }

        //! Closes the currently open file
        void close(close_handler_type const& on_close)
        {
            if (!on_close.empty()) try
            {
                std::ostream strm(&m_File);
                on_close(boost::ref(strm));
            }
            catch (...)
            {
            }

            close_file();
        }

        //! Puts the buffered record to the file
        void put_record(std::string const& storage, open_handler_type const& on_open, close_handler_type const& on_close)
        {
            // Check if it's time to rotate the file
            if (m_File.is_open() &&
                (m_RotationSize < m_Written + static_cast< uintmax_t >(storage.size()) ||
                (m_RotationInterval > 0U && m_RotationInterval < static_cast< unsigned int >(std::time(NULL) - m_LastRotation))))
            {
                close(on_close);
            }

            if (!m_File.is_open())
            {
                open_file();
                if (!on_open.empty()) try
                {
                    std::ostream strm(&m_File);
                    on_open(boost::ref(strm));
                }
                catch (...)
                {
                }
            }

            // Now put the data to the file
            write_data(storage);
        }

        //! Writes data from storage to the file
        void write_data(std::string const& storage)
        {
            BOOST_STATIC_ASSERT(std::numeric_limits< std::streamsize >::is_specialized);
            const std::streamsize max_block_size = (std::numeric_limits< std::streamsize >::max)();
            register std::string::size_type size = storage.size();
            register std::string::const_pointer const end = storage.data() + size;
            while (size > 0)
            {
                register const std::streamsize block_size =
                    size > static_cast< std::string::size_type >(max_block_size) ? max_block_size : static_cast< std::streamsize >(size);
                register const std::streamsize written = m_File.sputn(end - size, block_size);
                if (written == 0)
                    break;
                size -= static_cast< std::string::size_type >(written);
                m_Written += static_cast< uintmax_t >(written);
            }
        }

    private:
        //  Copying and assignment are prohibited
        file_controller_base(file_controller_base const&);
        file_controller_base& operator= (file_controller_base const&);

        //! The function opens a new file
        virtual void open_file() = 0;
        //! The function closes the file
        virtual void close_file() = 0;
    };

    //! The file controller implementation
    template< typename CharT >
    class file_controller :
        public file_controller_base
    {
        //! Base type
        typedef file_controller_base base_type;

    public:
        typedef CharT char_type;
        typedef std::basic_string< char_type > string_type;

#ifndef BOOST_FILESYSTEM_NARROW_ONLY
        typedef typename mpl::if_<
            is_same< char_type, char >,
            filesystem::path,
            filesystem::wpath
        >::type path_type;
#else
        typedef filesystem::path path_type;
#endif // BOOST_FILESYSTEM_NARROW_ONLY

    private:
        //! Date and time formatter
        class date_and_time_formatter
        {
            typedef date_time::time_facet< posix_time::ptime, char_type > time_facet_type;
            time_facet_type* m_pFacet;
            std::basic_ostringstream< char_type > m_Stream;
            const string_type m_EmptyString;

        public:
            //! Constructor
            date_and_time_formatter() : m_pFacet(NULL)
            {
                std::auto_ptr< time_facet_type > pFacet(new time_facet_type());
                m_pFacet = pFacet.get();
                std::locale loc(m_Stream.getloc(), m_pFacet);
                pFacet.release();
                m_Stream.imbue(loc);
            }
            //! The method formats the current date and time according to the format string str and writes the result into it
            void format(string_type& str)
            {
                m_pFacet->format(str.c_str());
                m_Stream.str(m_EmptyString);
                m_Stream << boost::log::attributes::local_time_traits::get_clock();
                if (m_Stream.good())
                    str = m_Stream.str();
                else
                    m_Stream.clear();
            }
        };

    private:
        //! The directory path
        path_type m_Path;
        //! The file name pattern
        string_type m_Pattern;

        //! The position in the pattern where the file counter placeholder is
        typename string_type::size_type m_FileCounterPosition;
        //! The file counter format
        string_type m_FileCounterFormat;

        //! The file counter
        unsigned int m_FileCounter;

        //! The date and time stamp formatter
        optional< date_and_time_formatter > m_DateFormatter;

    public:
        //! Constructor
        explicit file_controller(
            path_type const& pattern,
            std::ios_base::openmode mode,
            uintmax_t rot_size,
            unsigned int rot_int
        ) :
            base_type(mode, rot_size, rot_int),
            m_Path(pattern.branch_path()),
            m_Pattern(pattern.leaf()),
            m_FileCounterPosition(0),
            m_FileCounter(0)
        {
            typedef file_controller_traits< char_type > traits_t;
            typename string_type::const_iterator end = m_Pattern.end();
            typename string_type::const_iterator it = std::find(
                static_cast< typename string_type::const_iterator >(m_Pattern.begin()),
                end,
                static_cast< char_type >(traits_t::percent));
            unsigned int placeholder_count = 0;
            while (it != end)
            {
                ++placeholder_count;
                typename string_type::const_iterator placeholder_begin = it++;
                if (it == end)
                    break;
                typename string_type::const_iterator placeholder_end = parse_format_flag(it, end);

                if (placeholder_end != end)
                {
                    // We've found the file counter placeholder in the pattern
                    m_FileCounterFormat.assign(placeholder_begin, placeholder_end);
                    m_FileCounterFormat.push_back(
                        static_cast< typename string_type::value_type >(traits_t::digit_placeholder));
                    m_FileCounterPosition = placeholder_begin - m_Pattern.begin();
                    m_Pattern.erase(m_FileCounterPosition, placeholder_end - placeholder_begin + 1);
                    --placeholder_count;
                    break;
                }
            }

            // Yes, this is not totally correct, because the pattern can contain '%' which are not placeholders.
            // It just optimizes a bit in the obvious cases, so this is not crucial.
            if (placeholder_count > 0)
                m_DateFormatter = boost::in_place();
        }

        //! The function opens the file
        void open_file()
        {
            typedef file_controller_traits< char_type > traits_t;

            // Let's construct the new file name
            string_type file_name = m_Pattern;
            if (!m_FileCounterFormat.empty())
            {
                // Calculate the needed buffer length to print the counter.
                // The std::numeric_limits< unsigned int >::digits10 shows the maximum number of decimals
                // in the printed counter minus one.
                // Also, it is possible to specify width of the printed value, which may be 9 at maximum.
                // We should also count for a possible sign character (who would ever need it, however?) and
                // the terminating null character.
                enum { unsigned_int_length = std::numeric_limits< unsigned int >::digits10 + 1 };
                char_type counter_buf[(unsigned_int_length < 9 ? 9 : unsigned_int_length) + 2];
                int printed = traits_t::get_lazy_sprintf()(
                    counter_buf, sizeof(counter_buf) / sizeof(*counter_buf), m_FileCounterFormat.c_str(), m_FileCounter);
                if (printed > 0)
                {
                    file_name.insert(m_FileCounterPosition, counter_buf, printed);
                }
            }
            ++m_FileCounter;

            if (!!m_DateFormatter)
                m_DateFormatter->format(file_name);

            path_type full_file_name = m_Path / file_name;
            this->m_File.open(full_file_name, this->m_OpenMode);
            if (this->m_File.is_open())
                this->m_LastRotation = std::time(NULL);
            else
                boost::log::aux::throw_exception(std::runtime_error("failed to open file"));
        }

        //! The function closes the file
        void close_file()
        {
            if (this->m_File.is_open())
            {
                this->m_Written = 0;
                this->m_File.close();
            }
        }

    private:
        static typename string_type::const_iterator parse_format_flag(
            typename string_type::const_iterator& it, typename string_type::const_iterator end)
        {
            typedef file_controller_traits< char_type > traits_t;
            switch (*it)
            {
                case traits_t::plus:
                case traits_t::minus:
                case traits_t::space:
                case traits_t::zero:
                    // Format flag detected
                    if (++it == end)
                        return end;
                    if (traits_t::is_digit(*it))
                    {
                        // Format width detected
                        if (++it == end)
                            return end;
                    }
                default:
                    return parse_format_precision(it, end);
            }
        }
        static typename string_type::const_iterator parse_format_precision(
            typename string_type::const_iterator& it, typename string_type::const_iterator end)
        {
            typedef file_controller_traits< char_type > traits_t;
            switch (*it)
            {
                case traits_t::dot:
                    // Format precision detected
                    if (++it == end)
                        return end;
                    if (traits_t::is_digit(*it))
                    {
                        if (++it == end)
                            return end;
                    }
                    else
                        return end;

                    if (*it != traits_t::number_placeholder)
                        return end;

                case traits_t::number_placeholder:
                    // Placeholder detected
                    return it;

                default:
                    return end;
            }
        }
    };

} // namespace aux


/*!
 * \brief A rotating file output stream class
 *
 * The \c basic_rotating_ofstream class provides interface similar to standard
 * output stream except that it derives from the \c record_writer interface and requires
 * caller to invoke methods of the interface appropriately.
 *
 * The stream automatically manages underlying files depending on rotation criteria. The
 * stream also supports calling custom functional objects on closing and opening files.
 */
template<
    typename CharT,
    typename TraitsT = std::char_traits< CharT >,
    typename AllocatorT = std::allocator< CharT >
>
class basic_rotating_ofstream :
    public record_writer,
    public std::basic_ostream< CharT, TraitsT >
{
    //! Base class of the stream
    typedef std::basic_ostream< CharT, TraitsT > stream_base;

public:
    //! String type
    typedef std::basic_string< CharT, TraitsT, AllocatorT > string_type;

    //! The type of the handler of the event of opening of a new file
    typedef aux::file_controller_base::open_handler_type open_handler_type;
    //! The type of the handler of the event of closing of a file
    typedef aux::file_controller_base::close_handler_type close_handler_type;

private:

#ifndef BOOST_LOG_DOXYGEN_PASS

    //! Composite stream buffer type, that performs encoding translation and passes converted data to the file
    class composite_streambuf :
        private base_from_member< std::string >,
        public boost::log::aux::make_narrowing_ostringstreambuf<
            CharT,
            TraitsT,
            AllocatorT
        >::type
    {
    private:
        //! Narrowing buffer type
        typedef typename boost::log::aux::make_narrowing_ostringstreambuf<
            CharT,
            TraitsT,
            AllocatorT
        >::type narrowing_streambuf;

    private:
        //! The pointer to the file controller
        std::auto_ptr< aux::file_controller_base > m_pFileCtl;

    public:
        composite_streambuf() :
            narrowing_streambuf(base_from_member< std::string >::member)
        {
        }

        //! The function creates a file controller
        void open(filesystem::path const& pattern, std::ios_base::openmode mode, uintmax_t rot_size, unsigned int rot_int)
        {
            m_pFileCtl.reset(new aux::file_controller< char >(pattern, mode, rot_size, rot_int));
        }

#ifndef BOOST_FILESYSTEM_NARROW_ONLY
        //! The function creates a file controller
        void open(filesystem::wpath const& pattern, std::ios_base::openmode mode, uintmax_t rot_size, unsigned int rot_int)
        {
            m_pFileCtl.reset(new aux::file_controller< wchar_t >(pattern, mode, rot_size, rot_int));
        }
#endif // BOOST_FILESYSTEM_NARROW_ONLY

        //! The method closes the file
        void close(close_handler_type const& on_close)
        {
            boost::log::aux::cleanup_guard< std::string > _(storage());
            narrowing_streambuf::sync();

            if (m_pFileCtl.get())
            {
                m_pFileCtl->write_data(storage());
                m_pFileCtl->close(on_close);
                m_pFileCtl.reset();
            }
        }

        //! The method is called after all data of the record is written to the stream
        void on_end_record(open_handler_type const& on_open, close_handler_type const& on_close)
        {
            boost::log::aux::cleanup_guard< std::string > _(storage());
            narrowing_streambuf::sync();

            if (m_pFileCtl.get()) try
            {
                m_pFileCtl->put_record(storage(), on_open, on_close);
            }
            catch (...)
            {
            }
        }

    protected:
        //! Buffer synchronization implementation (flushes the attached file)
        virtual int sync()
        {
            // We actually flush file buffers, not the string buffers
            if (m_pFileCtl.get())
                return m_pFileCtl->flush();
            else
                return 0;
        }

    private:
        //  Copying and assignment are prohibited
        composite_streambuf(composite_streambuf const&);
        composite_streambuf& operator= (composite_streambuf const&);

        //! The function returns a reference to the storage string
        std::string& storage() { return base_from_member< std::string >::member; }
    };

#endif // BOOST_LOG_DOXYGEN_PASS

private:
    //! The stream buffer converts record data and puts it to the file
    composite_streambuf m_Buf;

    //! The handler in called when a new file is opened during rotation
    open_handler_type m_OpenHandler;
    //! The handler in called when a file is closed during rotation
    close_handler_type m_CloseHandler;

public:
    /*!
     * Default constructor. Creates the stream object in closed state
     * (i.e. no writing is possible until \c open is called).
     */
    basic_rotating_ofstream() : stream_base(NULL)
    {
        stream_base::init(&m_Buf);
    }

#ifndef BOOST_LOG_DOXYGEN_PASS

#ifndef BOOST_FILESYSTEM_NARROW_ONLY

#define BOOST_LOG_ROTATING_OFSTREAM_CTOR(z, it, data)\
    template< BOOST_PP_ENUM_PARAMS(it, typename T) >\
    basic_rotating_ofstream(filesystem::path const& pattern, BOOST_PP_ENUM_BINARY_PARAMS(it, T, const& arg)) :\
        stream_base(NULL)\
    {\
        stream_base::init(&m_Buf);\
        open_internal(pattern, (BOOST_PP_ENUM_PARAMS(it, arg)));\
    }\
    template< BOOST_PP_ENUM_PARAMS(it, typename T) >\
    basic_rotating_ofstream(filesystem::wpath const& pattern, BOOST_PP_ENUM_BINARY_PARAMS(it, T, const& arg)) :\
        stream_base(NULL)\
    {\
        stream_base::init(&m_Buf);\
        open_internal(pattern, (BOOST_PP_ENUM_PARAMS(it, arg)));\
    }\
    template< BOOST_PP_ENUM_PARAMS(it, typename T) >\
    void open(filesystem::path const& pattern, BOOST_PP_ENUM_BINARY_PARAMS(it, T, const& arg))\
    {\
        close();\
        open_internal(pattern, (BOOST_PP_ENUM_PARAMS(it, arg)));\
    }\
    template< BOOST_PP_ENUM_PARAMS(it, typename T) >\
    void open(filesystem::wpath const& pattern, BOOST_PP_ENUM_BINARY_PARAMS(it, T, const& arg))\
    {\
        close();\
        open_internal(pattern, (BOOST_PP_ENUM_PARAMS(it, arg)));\
    }

#else

#define BOOST_LOG_ROTATING_OFSTREAM_CTOR(z, it, data)\
    template< BOOST_PP_ENUM_PARAMS(it, typename T) >\
    basic_rotating_ofstreambuf(filesystem::path const& pattern, BOOST_PP_ENUM_BINARY_PARAMS(it, T, const& arg)) :\
        stream_base(NULL)\
    {\
        stream_base::init(&m_Buf);\
        open_internal(pattern, (BOOST_PP_ENUM_PARAMS(it, arg)));\
    }\
    template< BOOST_PP_ENUM_PARAMS(it, typename T) >\
    void open(filesystem::path const& pattern, BOOST_PP_ENUM_BINARY_PARAMS(it, T, const& arg))\
    {\
        close();\
        open_internal(pattern, (BOOST_PP_ENUM_PARAMS(it, arg)));\
    }

#endif // BOOST_FILESYSTEM_NARROW_ONLY

    BOOST_PP_REPEAT_FROM_TO(1, 4, BOOST_LOG_ROTATING_OFSTREAM_CTOR, ~)

#undef BOOST_LOG_ROTATING_OFSTREAM_CTOR

#else // BOOST_LOG_DOXYGEN_PASS

    /*!
     * Constructor. Creates the stream in open state.
     *
     * \note The construction does not lead to immediate opening or creating a file. The file will be opened
     *       on the first actual writing operation.
     *
     * \param pattern File name pattern. May contain <tt>%%N</tt> placeholder, which will be replaced with the
     *                file counter, or any of Boost.DateTime formatting placeholders for date and/or time.
     * \param args A set of named arguments. The following arguments are supported:
     *             \li \c rotation_size - file size to change the file at
     *             \li \c rotation_interval - time interval at which to rotate the file
     *             \li \c open_mode - a set of flags that describe the way to open files. Should be composed of
     *                    <tt>std::ios_base::openmode</tt> flags.
     */
    template< typename... ArgsT >
    basic_rotating_ofstreambuf(filesystem::path const& pattern, ArgsT... const& args);
    /*!
     * Constructor. Creates the stream in open state.
     *
     * \note The construction does not lead to immediate opening or creating a file. The file will be opened
     *       on the first actual writing operation.
     *
     * \param pattern File name pattern. May contain <tt>%%N</tt> placeholder, which will be replaced with the
     *                file counter, or any of Boost.DateTime formatting placeholders for date and/or time.
     * \param args A set of named arguments. The following arguments are supported:
     *             \li \c rotation_size - file size to change the file at
     *             \li \c rotation_interval - time interval at which to rotate the file
     *             \li \c open_mode - a set of flags that describe the way to open files. Should be composed of
     *                    <tt>std::ios_base::openmode</tt> flags.
     */
    template< typename... ArgsT >
    basic_rotating_ofstreambuf(filesystem::wpath const& pattern, ArgsT... const& args);
    /*!
     * Moves the stream into the open state.
     *
     * \note The method call does not lead to immediate opening or creating a file. The file will be opened
     *       on the first actual writing operation.
     *
     * \param pattern File name pattern. May contain <tt>%%N</tt> placeholder, which will be replaced with the
     *                file counter, or any of Boost.DateTime formatting placeholders for date and/or time.
     * \param args A set of named arguments. The following arguments are supported:
     *             \li \c rotation_size - file size to change the file at
     *             \li \c rotation_interval - time interval at which to rotate the file
     *             \li \c open_mode - a set of flags that describe the way to open files. Should be composed of
     *                    <tt>std::ios_base::openmode</tt> flags.
     */
    template< typename... ArgsT >
    void open(filesystem::path const& pattern, ArgsT... const& args);
    /*!
     * Moves the stream into the open state.
     *
     * \note The method call does not lead to immediate opening or creating a file. The file will be opened
     *       on the first actual writing operation.
     *
     * \param pattern File name pattern. May contain <tt>%%N</tt> placeholder, which will be replaced with the
     *                file counter, or any of Boost.DateTime formatting placeholders for date and/or time.
     * \param args A set of named arguments. The following arguments are supported:
     *             \li \c rotation_size - file size to change the file at
     *             \li \c rotation_interval - time interval at which to rotate the file
     *             \li \c open_mode - a set of flags that describe the way to open files. Should be composed of
     *                    <tt>std::ios_base::openmode</tt> flags.
     */
    template< typename... ArgsT >
    void open(filesystem::wpath const& pattern, ArgsT... const& args);

#endif // BOOST_LOG_DOXYGEN_PASS

    /*!
     * Destructor. Automatically closes the file, if it is open.
     */
    ~basic_rotating_ofstream()
    {
        try
        {
            close();
        }
        catch (...)
        {
        }
    }

    /*!
     * The method closes the file, if open, and moves the stream into the closed state.
     */
    void close()
    {
        m_Buf.close(m_CloseHandler);
    }

    virtual void on_end_record()
    {
        m_Buf.on_end_record(m_OpenHandler, m_CloseHandler);
    }

    /*!
     * Sets a new handler for opening a new file
     *
     * \param handler Functional object to receive notifications
     */
    void set_open_handler(open_handler_type const& handler)
    {
        m_OpenHandler = handler;
    }
    /*!
     * Removes the handler for opening a new file
     */
    void clear_open_handler()
    {
        m_OpenHandler.clear();
    }
    /*!
     * Sets a new handler for closing a file
     *
     * \param handler Functional object to receive notifications
     */
    void set_close_handler(close_handler_type const& handler)
    {
        m_CloseHandler = handler;
    }
    /*!
     * Removes the handler for closing a file
     */
    void clear_close_handler()
    {
        m_CloseHandler.clear();
    }

private:
#ifndef BOOST_LOG_DOXYGEN_PASS
    //  Copying and assignment are closed
    basic_rotating_ofstream(basic_rotating_ofstream const&);
    basic_rotating_ofstream& operator= (basic_rotating_ofstream const&);

    //! An internal function to handle construction args
    template< typename PathT, typename ArgsT >
    void open_internal(PathT const& pattern, ArgsT const& args)
    {
        m_Buf.open(
            pattern,
            args[keywords::open_mode | (std::ios_base::out | std::ios_base::trunc)],
            args[keywords::rotation_size | ~static_cast< uintmax_t >(0)],
            args[keywords::rotation_interval | static_cast< unsigned int >(0)]);
    }
#endif // BOOST_LOG_DOXYGEN_PASS
};

#ifdef BOOST_LOG_USE_CHAR
typedef basic_rotating_ofstream< char > rotating_ofstream;      //!< Convenience typedef for narrow-character logging
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
typedef basic_rotating_ofstream< wchar_t > rotating_wofstream;  //!< Convenience typedef for wide-character logging
#endif

} // namespace log

} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif // _MSC_VER

#endif // BOOST_LOG_SINKS_ROTATING_OFSTREAM_HPP_INCLUDED_
