/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   syslog_backend.cpp
 * \author Andrey Semashev
 * \date   08.01.2008
 *
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#include "windows_version.hpp"
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/socket_base.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/asio/ip/udp.hpp>
#include <boost/asio/ip/address.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/system/error_code.hpp>
#include <boost/date_time/c_time.hpp>
#include <boost/compatibility/cpp_c_headers/ctime>
#include <boost/log/sinks/syslog_backend.hpp>
#include <boost/log/detail/singleton.hpp>
#include <boost/log/detail/snprintf.hpp>
#include <boost/log/detail/throw_exception.hpp>
#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#endif

#ifdef BOOST_LOG_USE_NATIVE_SYSLOG
#include <syslog.h>
#endif // BOOST_LOG_USE_NATIVE_SYSLOG

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace sinks {

namespace syslog {

    //  Syslog record levels
    BOOST_LOG_EXPORT extern const level_t emergency = { 0 };
    BOOST_LOG_EXPORT extern const level_t alert = { 1 };
    BOOST_LOG_EXPORT extern const level_t critical = { 2 };
    BOOST_LOG_EXPORT extern const level_t error = { 3 };
    BOOST_LOG_EXPORT extern const level_t warning = { 4 };
    BOOST_LOG_EXPORT extern const level_t notice = { 5 };
    BOOST_LOG_EXPORT extern const level_t info = { 6 };
    BOOST_LOG_EXPORT extern const level_t debug = { 7 };

    //  Syslog facility codes
    BOOST_LOG_EXPORT extern const facility_t kernel = { 0 * 8 };
    BOOST_LOG_EXPORT extern const facility_t user = { 1 * 8 };
    BOOST_LOG_EXPORT extern const facility_t mail = { 2 * 8 };
    BOOST_LOG_EXPORT extern const facility_t daemon = { 3 * 8 };
    BOOST_LOG_EXPORT extern const facility_t security0 = { 4 * 8 };
    BOOST_LOG_EXPORT extern const facility_t syslogd = { 5 * 8 };
    BOOST_LOG_EXPORT extern const facility_t printer = { 6 * 8 };
    BOOST_LOG_EXPORT extern const facility_t news = { 7 * 8 };
    BOOST_LOG_EXPORT extern const facility_t uucp = { 8 * 8 };
    BOOST_LOG_EXPORT extern const facility_t clock0 = { 9 * 8 };
    BOOST_LOG_EXPORT extern const facility_t security1 = { 10 * 8 };
    BOOST_LOG_EXPORT extern const facility_t ftp = { 11 * 8 };
    BOOST_LOG_EXPORT extern const facility_t ntp = { 12 * 8 };
    BOOST_LOG_EXPORT extern const facility_t log_audit = { 13 * 8 };
    BOOST_LOG_EXPORT extern const facility_t log_alert = { 14 * 8 };
    BOOST_LOG_EXPORT extern const facility_t clock1 = { 15 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local0 = { 16 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local1 = { 17 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local2 = { 18 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local3 = { 19 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local4 = { 20 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local5 = { 21 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local6 = { 22 * 8 };
    BOOST_LOG_EXPORT extern const facility_t local7 = { 23 * 8 };

} // namespace syslog

////////////////////////////////////////////////////////////////////////////////
//! Syslog sink backend implementation
////////////////////////////////////////////////////////////////////////////////
template< typename CharT >
struct basic_syslog_backend< CharT >::implementation
{
#ifdef BOOST_LOG_USE_NATIVE_SYSLOG
    struct native;
#endif // BOOST_LOG_USE_NATIVE_SYSLOG
    struct udp_socket_based;

    //! Level mapper
    severity_mapper_type m_LevelMapper;

    //! Logging facility (portable or native, depending on the backend implementation)
    const int m_Facility;

    //! Constructor
    explicit implementation(int facility) :
        m_Facility(facility)
    {
    }
    //! Virtual destructor
    virtual ~implementation() {}

    //! The method sends the formatted message to the syslog host
    virtual void send(syslog::level_t level, target_string_type const& formatted_message) = 0;
};


////////////////////////////////////////////////////////////////////////////////
//  Native syslog API support
////////////////////////////////////////////////////////////////////////////////

#ifdef BOOST_LOG_USE_NATIVE_SYSLOG

namespace {

    //! Syslog service initializer (implemented as a weak singleton)
#if !defined(BOOST_LOG_NO_THREADS)
    class native_syslog_initializer :
        private log::aux::lazy_singleton< native_syslog_initializer, mutex >
#else
    class native_syslog_initializer
#endif
    {
#if !defined(BOOST_LOG_NO_THREADS)
        friend class log::aux::lazy_singleton< native_syslog_initializer, mutex >;
        typedef log::aux::lazy_singleton< native_syslog_initializer, mutex > mutex_holder;
#endif

    public:
        native_syslog_initializer()
        {
            ::openlog("", 0, LOG_USER);
        }
        ~native_syslog_initializer()
        {
            ::closelog();
        }

        static shared_ptr< native_syslog_initializer > get_instance()
        {
#if !defined(BOOST_LOG_NO_THREADS)
            lock_guard< mutex > _(mutex_holder::get());
#endif
            static weak_ptr< native_syslog_initializer > instance;
            shared_ptr< native_syslog_initializer > p(instance.lock());
            if (!p)
            {
                p = boost::make_shared< native_syslog_initializer >();
                instance = p;
            }
            return p;
        }
    };

} // namespace

template< typename CharT >
struct basic_syslog_backend< CharT >::implementation::native :
    public implementation
{
    //! Reference to the syslog service initializer
    const shared_ptr< native_syslog_initializer > m_pSyslogInitializer;

    //! Constructor
    explicit native(syslog::facility_t const& facility) :
        implementation(convert_facility(facility)),
        m_pSyslogInitializer(native_syslog_initializer::get_instance())
    {
    }

    //! The method sends the formatted message to the syslog host
    void send(syslog::level_t level, target_string_type const& formatted_message)
    {
        int native_level;
        switch (static_cast< int >(level))
        {
        case 0:
            native_level = LOG_EMERG; break;
        case 1:
            native_level = LOG_ALERT; break;
        case 2:
            native_level = LOG_CRIT; break;
        case 3:
            native_level = LOG_ERR; break;
        case 4:
            native_level = LOG_WARNING; break;
        case 5:
            native_level = LOG_NOTICE; break;
        case 7:
            native_level = LOG_DEBUG; break;
        default:
            native_level = LOG_INFO; break;
        }

        ::syslog(LOG_MAKEPRI(this->m_Facility, native_level), "%s", formatted_message.c_str());
    }

private:
    //! The function converts portable facility codes to the native codes
    static int convert_facility(syslog::facility_t const& facility)
    {
        // POSIX does not specify anything except for LOG_USER and LOG_LOCAL*
        #ifndef LOG_KERN
        #define LOG_KERN LOG_USER
        #endif
        #ifndef LOG_DAEMON
        #define LOG_DAEMON LOG_KERN
        #endif
        #ifndef LOG_MAIL
        #define LOG_MAIL LOG_USER
        #endif
        #ifndef LOG_AUTH
        #define LOG_AUTH LOG_DAEMON
        #endif
        #ifndef LOG_SYSLOG
        #define LOG_SYSLOG LOG_DAEMON
        #endif
        #ifndef LOG_LPR
        #define LOG_LPR LOG_DAEMON
        #endif
        #ifndef LOG_NEWS
        #define LOG_NEWS LOG_USER
        #endif
        #ifndef LOG_UUCP
        #define LOG_UUCP LOG_USER
        #endif
        #ifndef LOG_CRON
        #define LOG_CRON LOG_DAEMON
        #endif
        #ifndef LOG_AUTHPRIV
        #define LOG_AUTHPRIV LOG_AUTH
        #endif
        #ifndef LOG_FTP
        #define LOG_FTP LOG_DAEMON
        #endif

        static const int native_facilities[24] =
        {
            LOG_KERN,
            LOG_USER,
            LOG_MAIL,
            LOG_DAEMON,
            LOG_AUTH,
            LOG_SYSLOG,
            LOG_LPR,
            LOG_NEWS,
            LOG_UUCP,
            LOG_CRON,
            LOG_AUTHPRIV,
            LOG_FTP,

            // reserved values
            LOG_USER,
            LOG_USER,
            LOG_USER,
            LOG_USER,

            LOG_LOCAL0,
            LOG_LOCAL1,
            LOG_LOCAL2,
            LOG_LOCAL3,
            LOG_LOCAL4,
            LOG_LOCAL5,
            LOG_LOCAL6,
            LOG_LOCAL7
        };

        register std::size_t n = static_cast< int >(facility) / 8;
        BOOST_ASSERT(n < sizeof(native_facilities) / sizeof(*native_facilities));
        return native_facilities[n];
    }
};

#endif // BOOST_LOG_USE_NATIVE_SYSLOG


////////////////////////////////////////////////////////////////////////////////
//  Socket-based implementation
////////////////////////////////////////////////////////////////////////////////
namespace {

    //! The shared UDP socket
    struct syslog_udp_socket
    {
    private:
        //! The socket primitive
        asio::ip::udp::socket m_Socket;

    public:
        //! The constructor creates a socket bound to the specified local address and port
        explicit syslog_udp_socket(asio::io_service& service, asio::ip::udp const& protocol, asio::ip::udp::endpoint const& local_address) :
            m_Socket(service)
        {
            m_Socket.open(protocol);
            m_Socket.set_option(asio::socket_base::reuse_address(true));
            m_Socket.bind(local_address);
        }
        //! The destructor closes the socket
        ~syslog_udp_socket()
        {
            boost::system::error_code ec;
            m_Socket.shutdown(asio::socket_base::shutdown_both, ec);
            m_Socket.close(ec);
        }

        //! The method sends the syslog message to the specified endpoint
        void send_message(int pri, const char* local_host_name, asio::ip::udp::endpoint const& target, const char* message);

    private:
        syslog_udp_socket(syslog_udp_socket const&);
        syslog_udp_socket& operator= (syslog_udp_socket const&);
    };

    //! The class contains the UDP service for syslog sockets to function
    class syslog_udp_service :
        public log::aux::lazy_singleton< syslog_udp_service, shared_ptr< syslog_udp_service > >
    {
        friend class log::aux::lazy_singleton< syslog_udp_service, shared_ptr< syslog_udp_service > >;
        typedef log::aux::lazy_singleton< syslog_udp_service, shared_ptr< syslog_udp_service > > base_type;

    public:
        //! The core IO service instance
        asio::io_service m_IOService;
        //! The local host name to put into log message
        std::string m_LocalHostName;

#if !defined(BOOST_LOG_NO_THREADS)
        //! A synchronization primitive to protect the host name resolver
        mutex m_Mutex;
#endif // !defined(BOOST_LOG_NO_THREADS)
        //! The resolver is used to acquire connection endpoints
        asio::ip::udp::resolver m_HostNameResolver;

    private:
        //! Default constructor
        syslog_udp_service() :
            m_HostNameResolver(m_IOService)
        {
            boost::system::error_code err;
            m_LocalHostName = asio::ip::host_name(err);
        }
        //! Initializes the singleton instance
        static void init_instance()
        {
            base_type::get_instance().reset(new syslog_udp_service());
        }
    };

    //! The method sends the syslog message to the specified endpoint
    void syslog_udp_socket::send_message(
        int pri, const char* local_host_name, asio::ip::udp::endpoint const& target, const char* message)
    {
        std::time_t t = std::time(NULL);
        std::tm ts;
        std::tm* time_stamp = boost::date_time::c_time::localtime(&t, &ts);

        // Month will have to be injected separately, as involving locale won't do here
        static const char months[12][4] =
        {
            "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
        };

        // The packet size is mandated in RFC3164, plus one for the terminating zero
        char packet[1025];
        std::size_t packet_size = boost::log::aux::snprintf(
            packet,
            sizeof(packet),
            "<%d> %s % 2d %02d:%02d:%02d %s %s",
            pri,
            months[time_stamp->tm_mon],
            time_stamp->tm_mday,
            time_stamp->tm_hour,
            time_stamp->tm_min,
            time_stamp->tm_sec,
            local_host_name,
            message
        );

        m_Socket.send_to(asio::buffer(packet, packet_size), target);
    }

} // namespace

template< typename CharT >
struct basic_syslog_backend< CharT >::implementation::udp_socket_based :
    public implementation
{
    //! Protocol to be used
    asio::ip::udp m_Protocol;
    //! Pointer to the list of sockets
    shared_ptr< syslog_udp_service > m_pService;
    //! Pointer to the socket being used
    std::auto_ptr< syslog_udp_socket > m_pSocket;
    //! The target host to send packets to
    asio::ip::udp::endpoint m_TargetHost;

    //! Constructor
    explicit udp_socket_based(syslog::facility_t const& facility, asio::ip::udp const& protocol) :
        implementation(facility),
        m_Protocol(protocol),
        m_pService(syslog_udp_service::get())
    {
        if (m_Protocol == asio::ip::udp::v4())
        {
            m_TargetHost = asio::ip::udp::endpoint(asio::ip::address_v4(0x7F000001), 514); // 127.0.0.1:514
        }
        else
        {
            // ::1, port 514
            asio::ip::address_v6::bytes_type addr;
            std::fill_n(addr.c_array(), addr.size() - 1, static_cast< unsigned char >(0));
            addr[addr.size() - 1] = 1;
            m_TargetHost = asio::ip::udp::endpoint(asio::ip::address_v6(addr), 514);
        }
    }

    //! The method sends the formatted message to the syslog host
    void send(syslog::level_t level, target_string_type const& formatted_message)
    {
        if (!m_pSocket.get())
        {
            asio::ip::udp::endpoint any_local_address;
            m_pSocket.reset(new syslog_udp_socket(m_pService->m_IOService, m_Protocol, any_local_address));
        }

        m_pSocket->send_message(
            this->m_Facility | static_cast< int >(level),
            m_pService->m_LocalHostName.c_str(),
            m_TargetHost,
            formatted_message.c_str());
    }
};


////////////////////////////////////////////////////////////////////////////////
//  Sink backend implementation
////////////////////////////////////////////////////////////////////////////////
//! Destructor
template< typename CharT >
basic_syslog_backend< CharT >::~basic_syslog_backend()
{
    delete m_pImpl;
}

//! The method installs the function object that maps application severity levels to Syslog levels
template< typename CharT >
void basic_syslog_backend< CharT >::set_severity_mapper(severity_mapper_type const& mapper)
{
    m_pImpl->m_LevelMapper = mapper;
}

//! The method writes the message to the sink
template< typename CharT >
void basic_syslog_backend< CharT >::do_consume(
    record_type const& record, target_string_type const& formatted_message)
{
    m_pImpl->send(
        m_pImpl->m_LevelMapper.empty() ? syslog::info : m_pImpl->m_LevelMapper(record.attribute_values()),
        formatted_message);
}


//! The method creates the backend implementation
template< typename CharT >
typename basic_syslog_backend< CharT >::implementation*
basic_syslog_backend< CharT >::construct(
    syslog::facility_t facility, syslog::impl_types use_impl, ip_versions ip_version)
{
#ifdef BOOST_LOG_USE_NATIVE_SYSLOG
    if (use_impl == syslog::native)
    {
        typedef typename implementation::native native_impl;
        return new native_impl(facility);
    }
#endif // BOOST_LOG_USE_NATIVE_SYSLOG

    typedef typename implementation::udp_socket_based udp_socket_based_impl;
    switch (ip_version)
    {
    case v4:
        return new udp_socket_based_impl(facility, asio::ip::udp::v4());
    case v6:
        return new udp_socket_based_impl(facility, asio::ip::udp::v6());
    default:
        boost::log::aux::throw_exception(std::invalid_argument("incorrect IP version specified"));
    }
}

//! The method sets the local address which log records will be sent from.
template< typename CharT >
void basic_syslog_backend< CharT >::set_local_address(std::string const& addr, unsigned short port)
{
    typedef typename implementation::udp_socket_based udp_socket_based_impl;
    if (udp_socket_based_impl* impl = dynamic_cast< udp_socket_based_impl* >(m_pImpl))
    {
        char service_name[std::numeric_limits< int >::digits10 + 3];
        boost::log::aux::snprintf(service_name, sizeof(service_name), "%d", static_cast< int >(port));
        asio::ip::udp::resolver::query q(
            impl->m_Protocol,
            addr,
            service_name,
            asio::ip::resolver_query_base::address_configured | asio::ip::resolver_query_base::passive);
        asio::ip::udp::endpoint local_address;
        
        {
#if !defined(BOOST_LOG_NO_THREADS)
            lock_guard< mutex > _(impl->m_pService->m_Mutex);
#endif

            local_address = *impl->m_pService->m_HostNameResolver.resolve(q);
        }

        impl->m_pSocket.reset(new syslog_udp_socket(impl->m_pService->m_IOService, impl->m_Protocol, local_address));
    }
}
//! The method sets the local address which log records will be sent from.
template< typename CharT >
void basic_syslog_backend< CharT >::set_local_address(boost::asio::ip::address const& addr, unsigned short port)
{
    typedef typename implementation::udp_socket_based udp_socket_based_impl;
    if (udp_socket_based_impl* impl = dynamic_cast< udp_socket_based_impl* >(m_pImpl))
    {
        impl->m_pSocket.reset(new syslog_udp_socket(
            impl->m_pService->m_IOService, impl->m_Protocol, asio::ip::udp::endpoint(addr, port)));
    }
}

//! The method sets the address of the remote host where log records will be sent to.
template< typename CharT >
void basic_syslog_backend< CharT >::set_target_address(std::string const& addr, unsigned short port)
{
    typedef typename implementation::udp_socket_based udp_socket_based_impl;
    if (udp_socket_based_impl* impl = dynamic_cast< udp_socket_based_impl* >(m_pImpl))
    {
        char service_name[std::numeric_limits< int >::digits10 + 3];
        boost::log::aux::snprintf(service_name, sizeof(service_name), "%d", static_cast< int >(port));
        asio::ip::udp::resolver::query q(impl->m_Protocol, addr, service_name, asio::ip::resolver_query_base::address_configured);
        asio::ip::udp::endpoint remote_address;

        {
#if !defined(BOOST_LOG_NO_THREADS)
            lock_guard< mutex > _(impl->m_pService->m_Mutex);
#endif

            remote_address = *impl->m_pService->m_HostNameResolver.resolve(q);
        }

        impl->m_TargetHost = remote_address;
    }
}
//! The method sets the address of the remote host where log records will be sent to.
template< typename CharT >
void basic_syslog_backend< CharT >::set_target_address(boost::asio::ip::address const& addr, unsigned short port)
{
    typedef typename implementation::udp_socket_based udp_socket_based_impl;
    if (udp_socket_based_impl* impl = dynamic_cast< udp_socket_based_impl* >(m_pImpl))
    {
        impl->m_TargetHost = asio::ip::udp::endpoint(addr, port);
    }
}

//! Explicitly instantiate sink implementation
#ifdef BOOST_LOG_USE_CHAR
template class basic_syslog_backend< char >;
#endif
#ifdef BOOST_LOG_USE_WCHAR_T
template class basic_syslog_backend< wchar_t >;
#endif

} // namespace sinks

} // namespace log

} // namespace boost
