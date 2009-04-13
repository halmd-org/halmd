/*!
 * (C) 2007 Andrey Semashev
 *
 * Use, modification and distribution is subject to the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * \file   filter_parser.cpp
 * \author Andrey Semashev
 * \date   31.03.2008
 *
 * \brief  This header is the Boost.Log library implementation, see the library documentation
 *         at http://www.boost.org/libs/log/doc/log.html.
 */

#ifndef BOOST_LOG_NO_SETTINGS_PARSERS_SUPPORT

#include <map>
#include <stack>
#include <string>
#include <stdexcept>

#if !defined(BOOST_LOG_NO_THREADS) && !defined(BOOST_SPIRIT_THREADSAFE)
#define BOOST_SPIRIT_THREADSAFE
#endif // !defined(BOOST_LOG_NO_THREADS) && !defined(BOOST_SPIRIT_THREADSAFE)

#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/none.hpp>
#include <boost/optional.hpp>
#include <boost/utility/addressof.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_confix.hpp>
#include <boost/spirit/include/classic_escape_char.hpp>
#include <boost/log/core.hpp>
#include <boost/log/filters/basic_filters.hpp>
#include <boost/log/filters/attr.hpp>
#include <boost/log/filters/has_attr.hpp>
#include <boost/log/detail/singleton.hpp>
#include <boost/log/detail/functional.hpp>
#include <boost/log/detail/throw_exception.hpp>
#include <boost/log/utility/init/filter_parser.hpp>
#include <boost/log/utility/type_dispatch/standard_types.hpp>
#if !defined(BOOST_LOG_NO_THREADS)
#include <boost/log/detail/shared_lock_guard.hpp>
#include <boost/log/detail/light_rw_mutex.hpp>
#include <boost/thread/locks.hpp>
#endif // !defined(BOOST_LOG_NO_THREADS)
#include "parser_utils.hpp"

namespace boost {

namespace BOOST_LOG_NAMESPACE {

namespace {

//! The default filter factory that supports creating filters for the standard types (see utility/type_dispatch/standard_types.hpp)
template< typename CharT >
class default_filter_factory :
    public filter_factory< CharT >
{
    //! Base type
    typedef filter_factory< CharT > base_type;
    //! Self type
    typedef default_filter_factory< CharT > this_type;

    //  Type imports
    typedef typename base_type::char_type char_type;
    typedef typename base_type::string_type string_type;
    typedef typename base_type::filter_type filter_type;

    //! The callback for equality relation filter
    virtual filter_type on_equality_relation(string_type const& name, string_type const& arg)
    {
        return parse_argument< log::aux::equal_to >(name, arg);
    }
    //! The callback for inequality relation filter
    virtual filter_type on_inequality_relation(string_type const& name, string_type const& arg)
    {
        return parse_argument< log::aux::not_equal_to >(name, arg);
    }
    //! The callback for less relation filter
    virtual filter_type on_less_relation(string_type const& name, string_type const& arg)
    {
        return parse_argument< log::aux::less >(name, arg);
    }
    //! The callback for greater relation filter
    virtual filter_type on_greater_relation(string_type const& name, string_type const& arg)
    {
        return parse_argument< log::aux::greater >(name, arg);
    }
    //! The callback for less or equal relation filter
    virtual filter_type on_less_or_equal_relation(string_type const& name, string_type const& arg)
    {
        return parse_argument< log::aux::less_equal >(name, arg);
    }
    //! The callback for greater or equal relation filter
    virtual filter_type on_greater_or_equal_relation(string_type const& name, string_type const& arg)
    {
        return parse_argument< log::aux::greater_equal >(name, arg);
    }

    //! The callback for custom relation filter
    virtual filter_type on_custom_relation(string_type const& name, string_type const& rel, string_type const& arg)
    {
        typedef log::aux::char_constants< char_type > constants;
        if (rel == constants::begins_with_keyword())
            return filter_type(log::filters::attr< string_type >(name).begins_with(arg));
        else if (rel == constants::ends_with_keyword())
            return filter_type(log::filters::attr< string_type >(name).ends_with(arg));
        else if (rel == constants::contains_keyword())
            return filter_type(log::filters::attr< string_type >(name).contains(arg));
        else if (rel == constants::matches_keyword())
            return filter_type(log::filters::attr< string_type >(name).matches(arg));
        else
            boost::log::aux::throw_exception(std::runtime_error("the custom attribute relation is not supported"));

        // To get rid from compiler warnings
        BOOST_LOG_ASSUME(false);
        return filter_type();
    }

    //! The function parses the argument value for a binary relation and constructs the corresponding filter
    template< typename RelationT >
    static filter_type parse_argument(string_type const& name, string_type const& arg)
    {
        filter_type filter;
        spirit::classic::rule< spirit::classic::scanner< const char_type* > > r =
            spirit::classic::strict_real_p[bind(&this_type::BOOST_NESTED_TEMPLATE on_fp_argument< RelationT >, boost::cref(name), _1, boost::ref(filter))] |
            spirit::classic::int_p[bind(&this_type::BOOST_NESTED_TEMPLATE on_integral_argument< RelationT >, boost::cref(name), _1, boost::ref(filter))] |
            (+spirit::classic::print_p)[bind(&this_type::BOOST_NESTED_TEMPLATE on_string_argument< RelationT >, boost::cref(name), _1, _2, boost::ref(filter))];

        if (!spirit::classic::parse(arg.c_str(), r).full || filter.empty())
            boost::log::aux::throw_exception(std::runtime_error("failed to parse relation operand"));

        return filter;
    }

    template< typename RelationT >
    static void on_integral_argument(string_type const& name, long val, filter_type& filter)
    {
        filter = log::filters::attr<
            log::integral_types
        >(name).satisfies(log::aux::bind2nd(RelationT(), val));
    }
    template< typename RelationT >
    static void on_fp_argument(string_type const& name, double val, filter_type& filter)
    {
        filter = log::filters::attr<
            log::floating_point_types
        >(name).satisfies(log::aux::bind2nd(RelationT(), val));
    }
    template< typename RelationT >
    static void on_string_argument(string_type const& name, const char_type* b, const char_type* e, filter_type& filter)
    {
        filter = log::filters::attr<
            string_type
        >(name).satisfies(log::aux::bind2nd(RelationT(), string_type(b, e)));
    }
};

//! Filter factories repository
template< typename CharT >
struct filters_repository :
    public log::aux::lazy_singleton< filters_repository< CharT > >
{
    typedef CharT char_type;
    typedef log::aux::lazy_singleton< filters_repository< char_type > > base_type;
    typedef std::basic_string< char_type > string_type;
    typedef filter_factory< char_type > filter_factory_type;
    typedef std::map< string_type, shared_ptr< filter_factory_type > > factories_map;

#if !defined(BOOST_MSVC) || _MSC_VER > 1310
    friend class log::aux::lazy_singleton< filters_repository< char_type > >;
#else
    friend class base_type;
#endif

#if !defined(BOOST_LOG_NO_THREADS)
    //! Synchronization mutex
    mutable log::aux::light_rw_mutex m_Mutex;
#endif
    //! The map of filter factories
    factories_map m_Map;
    //! Default factory
    mutable default_filter_factory< char_type > m_DefaultFactory;

    //! The method returns the filter factory for the specified attribute name
    filter_factory_type& get_factory(string_type const& name) const
    {
        typename factories_map::const_iterator it = m_Map.find(name);
        if (it != m_Map.end())
            return *it->second;
        else
            return m_DefaultFactory;
    }

private:
    filters_repository() {}
};

//! Filter parsing grammar
template< typename CharT >
struct filter_grammar :
    public spirit::classic::grammar< filter_grammar< CharT > >
{
    typedef CharT char_type;
    typedef std::basic_string< char_type > string_type;
    typedef typename basic_core< char_type >::filter_type filter_type;
    typedef log::aux::char_constants< char_type > constants;
    typedef filter_grammar< char_type > filter_grammar_type;
    typedef filter_factory< char_type > filter_factory_type;

    typedef filter_type (filter_factory_type::*comparison_relation_handler_t)(string_type const&, string_type const&);

    template< typename ScannerT >
    struct definition;

    //! Parsed attribute name
    mutable optional< string_type > m_AttributeName;
    //! The second operand of a relation
    mutable optional< string_type > m_Operand;
    //! The custom relation string
    mutable string_type m_CustomRelation;

    //! Filter subexpressions as they are parsed
    mutable std::stack< filter_type > m_Subexpressions;

    //! Reference to the filter being constructed
    filter_type& m_Filter;

    //! Constructor
    explicit filter_grammar(filter_type& f) : m_Filter(f) {}

    //! The method finalizes filter construction by flushing its internal data that may not have been put into the filter
    void flush() const
    {
        if (!m_Subexpressions.empty())
            m_Filter.swap(m_Subexpressions.top());
    }

    //! The operand string handler
    void on_operand(const char_type* begin, const char_type* end) const
    {
        // An attribute name should have been parsed at this time
        if (!m_AttributeName)
            boost::log::aux::throw_exception(std::runtime_error("Invalid filter definition: operand is not expected"));

        m_Operand = boost::in_place(begin, end);
    }

    //! The quoted string handler
    void on_quoted_string(const char_type* begin, const char_type* end) const
    {
        // An attribute name should have been parsed at this time
        if (!m_AttributeName)
            boost::log::aux::throw_exception(std::runtime_error("Invalid filter definition: quoted string operand is not expected"));

        // Cut off the quotes
        string_type str(++begin, --end);

        // Translate escape sequences
        constants::translate_escape_sequences(str);
        m_Operand = str;
    }
    //! The attribute name handler
    void on_attribute(const char_type* begin, const char_type* end) const
    {
        // In case if previous subexpression consisted only
        // from attribute name, like in "%Attribute1% & %Attribute2% > 1"
        make_has_attr();

        // Cut off the '%'
        m_AttributeName = boost::in_place(++begin, --end);
    }

    //! The negation operation handler
    void on_negation(const char_type* begin, const char_type* end) const
    {
        make_has_attr();
        if (!m_Subexpressions.empty())
        {
            m_Subexpressions.top() = !log::filters::wrap(m_Subexpressions.top());
        }
        else
        {
            // This would happen if a filter consists of a single '!'
            boost::log::aux::throw_exception(std::runtime_error("Filter parsing error:"
                " a negation operator applied to nothingness"));
        }
    }
    //! The comparison relation handler
    void on_comparison_relation(const char_type* begin, const char_type* end, comparison_relation_handler_t method) const
    {
        if (!!m_AttributeName && !!m_Operand)
        {
            filters_repository< char_type > const& repo = filters_repository< char_type >::get();
            filter_factory_type& factory = repo.get_factory(m_AttributeName.get());

            m_Subexpressions.push((factory.*method)(m_AttributeName.get(), m_Operand.get()));

            m_AttributeName = none;
            m_Operand = none;
        }
        else
        {
            // This should never happen
            boost::log::aux::throw_exception(std::logic_error("Filter parser internal error:"
                " the attribute name or subexpression operand is not set while trying to construct a subexpression"));
        }
    }

    //! The method saves the relation word into an internal string
    void set_custom_relation(const char_type* begin, const char_type* end) const
    {
        m_CustomRelation.assign(begin, end);
    }

    //! The custom relation handler for string operands
    void on_custom_relation(const char_type* begin, const char_type* end) const
    {
        if (!!m_AttributeName && !!m_Operand && !m_CustomRelation.empty())
        {
            filters_repository< char_type > const& repo = filters_repository< char_type >::get();
            filter_factory_type& factory = repo.get_factory(m_AttributeName.get());

            m_Subexpressions.push(factory.on_custom_relation(m_AttributeName.get(), m_CustomRelation, m_Operand.get()));

            m_AttributeName = none;
            m_Operand = none;
            m_CustomRelation.clear();
        }
        else
        {
            // This should never happen
            boost::log::aux::throw_exception(std::logic_error("Filter parser internal error:"
                " the attribute name or subexpression operand is not set while trying to construct a subexpression"));
        }
    }

    //! The boolean operation handler
    template< template< typename, typename > class OperationT >
    void on_operation(const char_type* begin, const char_type* end) const
    {
        if (!m_Subexpressions.empty())
        {
            filter_type right = m_Subexpressions.top();
            m_Subexpressions.pop();
            if (!m_Subexpressions.empty())
            {
                filter_type const& left = m_Subexpressions.top();
                typedef log::filters::flt_wrap< char_type, filter_type > wrap_t;
                m_Subexpressions.top() = OperationT< wrap_t, wrap_t >(wrap_t(left), wrap_t(right));
                return;
            }
        }

        // This should never happen
        boost::log::aux::throw_exception(std::logic_error("Filter parser internal error:"
            " the subexpression is not set while trying to construct a filter"));
    }

    //! The function is called when a full expression have finished parsing
    void on_expression_finished(const char_type* begin, const char_type* end) const
    {
        make_has_attr();
    }

private:
    //  Assignment and copying are prohibited
    filter_grammar(filter_grammar const&);
    filter_grammar& operator= (filter_grammar const&);

    //! The function converts the parsed attribute name into a has_attr filter
    void make_has_attr() const
    {
        if (!!m_AttributeName)
        {
            filters_repository< char_type > const& repo = filters_repository< char_type >::get();
            filter_factory_type& factory = repo.get_factory(m_AttributeName.get());
            m_Subexpressions.push(factory.on_exists_test(m_AttributeName.get()));
            m_AttributeName = none;
        }
    }
};

//! Grammar definition
template< typename CharT >
template< typename ScannerT >
struct filter_grammar< CharT >::definition
{
    //! Boost.Spirit rule type
    typedef spirit::classic::rule< ScannerT > rule_type;

    //! A simple mem_fn-like wrapper (a workaround for MSVC 7.1)
    struct handler
    {
        typedef void result_type;
        typedef void (filter_grammar_type::*fun_type)(const char_type*, const char_type*) const;
        handler(fun_type pf) : m_fun(pf) {}

        void operator() (filter_grammar_type const* p, const char_type* b, const char_type* e) const
        {
            (p->*m_fun)(b, e);
        }

    private:
        fun_type m_fun;
    };

    //! A parser for an attribute name in a single relation
    rule_type attr_name;
    //! A parser for an operand in a single relation
    rule_type operand;
    //! A parser for a single relation that consists of two operands and an operation between them
    rule_type relation;
    //! A parser for a custom relation word
    rule_type custom_relation;
    //! A parser for a term, which can be a relation, an expression in parenthesis or a negation thereof
    rule_type term;
    //! A parser for the complete filter expression that consists of one or several terms with boolean operations between them
    rule_type expression;

    //! Constructor
    definition(filter_grammar_type const& gram)
    {
        const filter_grammar_type* g = boost::addressof(gram);

        // MSVC 7.1 goes wild for some reason if we try to use bind or mem_fn directly
        // on some of these functions. The simple wrapper helps the compiler to deduce types correctly.
        handler on_string = &filter_grammar_type::on_quoted_string,
            on_oper = &filter_grammar_type::on_operand,
            on_attr = &filter_grammar_type::on_attribute;

        attr_name = (*spirit::classic::space_p) >> (
            // An attribute name in form %name%
            spirit::classic::confix_p(constants::char_percent, *(spirit::classic::print_p - constants::char_percent), constants::char_percent)
                [bind(on_attr, g, _1, _2)]
        );

        operand = (*spirit::classic::space_p) >> (
            // A quoted string with C-style escape sequences support
            spirit::classic::confix_p(constants::char_quote, *spirit::classic::c_escape_ch_p, constants::char_quote)
                [bind(on_string, g, _1, _2)] |
            // A single word, enclosed with white spaces. It cannot contain parenthesis, since is is used by the filter parser.
            (+(spirit::classic::print_p - spirit::classic::space_p - constants::char_paren_bracket_left - constants::char_paren_bracket_right))
                [bind(on_oper, g, _1, _2)]
        );

        // Custom relation is a keyword that may contain either alphanumeric characters or an underscore
        custom_relation = (+(spirit::classic::alnum_p | spirit::classic::ch_p(constants::char_underline)))
            [bind(&filter_grammar_type::set_custom_relation, g, _1, _2)];

        relation = attr_name >> (*spirit::classic::space_p) || (
            (spirit::classic::str_p(constants::not_equal_keyword()) >> operand)
                [bind(&filter_grammar_type::on_comparison_relation, g, _1, _2, &filter_factory_type::on_inequality_relation)] |
            (spirit::classic::str_p(constants::greater_or_equal_keyword()) >> operand)
                [bind(&filter_grammar_type::on_comparison_relation, g, _1, _2, &filter_factory_type::on_greater_or_equal_relation)] |
            (spirit::classic::str_p(constants::less_or_equal_keyword()) >> operand)
                [bind(&filter_grammar_type::on_comparison_relation, g, _1, _2, &filter_factory_type::on_less_or_equal_relation)] |
            (constants::char_equal >> operand)
                [bind(&filter_grammar_type::on_comparison_relation, g, _1, _2, &filter_factory_type::on_equality_relation)] |
            (constants::char_greater >> operand)
                [bind(&filter_grammar_type::on_comparison_relation, g, _1, _2, &filter_factory_type::on_greater_relation)] |
            (constants::char_less >> operand)
                [bind(&filter_grammar_type::on_comparison_relation, g, _1, _2, &filter_factory_type::on_less_relation)] |
            (custom_relation >> operand)
                [bind(&filter_grammar_type::on_custom_relation, g, _1, _2)]
        );

        handler on_neg = &filter_grammar_type::on_negation;

        term = (*spirit::classic::space_p) >> (
            (constants::char_paren_bracket_left >> expression >> (*spirit::classic::space_p) >> constants::char_paren_bracket_right) |
            ((spirit::classic::str_p(constants::not_keyword()) | constants::char_exclamation) >> term)[bind(on_neg, g, _1, _2)] |
            relation
        );

        handler on_and = &filter_grammar_type::BOOST_NESTED_TEMPLATE on_operation< log::filters::flt_and >,
            on_or = &filter_grammar_type::BOOST_NESTED_TEMPLATE on_operation< log::filters::flt_or >,
            on_finished = &filter_grammar_type::on_expression_finished;

        expression = (term >> (*spirit::classic::space_p) >> *(
            ((spirit::classic::str_p(constants::and_keyword()) | constants::char_and) >> term)
                [bind(on_and, g, _1, _2)] |
            ((spirit::classic::str_p(constants::or_keyword()) | constants::char_or) >> term)
                [bind(on_or, g, _1, _2)]
            ) >> (*spirit::classic::space_p)
        )[bind(on_finished, g, _1, _2)];
    }

    //! Accessor for the filter rule
    rule_type const& start() const { return expression; }
};

} // namespace

//! The function registers a filter factory object for the specified attribute name
template< typename CharT >
void register_filter_factory(const CharT* attr_name, shared_ptr< filter_factory< CharT > > const& factory)
{
    std::basic_string< CharT > name(attr_name);
    filters_repository< CharT >& repo = filters_repository< CharT >::get();

#if !defined(BOOST_LOG_NO_THREADS)
    lock_guard< log::aux::light_rw_mutex > _(repo.m_Mutex);
#endif
    repo.m_Map[name] = factory;
}

//! The function parses a filter from the string
template< typename CharT >
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
typename basic_core< CharT >::filter_type
#else
function1< bool, basic_attribute_values_view< CharT > const& >
#endif
parse_filter(const CharT* begin, const CharT* end)
{
    typedef CharT char_type;
    typedef typename basic_core< char_type >::filter_type filter_type;

#if !defined(BOOST_LOG_NO_THREADS)
    filters_repository< CharT >& repo = filters_repository< CharT >::get();
    log::aux::shared_lock_guard< log::aux::light_rw_mutex > _(repo.m_Mutex);
#endif

    filter_type filt;
    filter_grammar< char_type > gram(filt);
    if (!spirit::classic::parse(begin, end, gram).full)
        boost::log::aux::throw_exception(std::runtime_error("Could not parse the filter"));
    gram.flush();

    return filt;
}

#ifdef BOOST_LOG_USE_CHAR

template BOOST_LOG_EXPORT
void register_filter_factory(const char* name, shared_ptr< filter_factory< char > > const& factory);

template BOOST_LOG_EXPORT
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
basic_core< char >::filter_type
#else
function1< bool, basic_attribute_values_view< char > const& >
#endif // BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
parse_filter< char >(const char* begin, const char* end);

#endif // BOOST_LOG_USE_CHAR

#ifdef BOOST_LOG_USE_WCHAR_T

template BOOST_LOG_EXPORT
void register_filter_factory(const wchar_t* name, shared_ptr< filter_factory< wchar_t > > const& factory);

template BOOST_LOG_EXPORT
#ifndef BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
basic_core< wchar_t >::filter_type
#else
function1< bool, basic_attribute_values_view< wchar_t > const& >
#endif // BOOST_LOG_BROKEN_TEMPLATE_DEFINITION_MATCHING
parse_filter< wchar_t >(const wchar_t* begin, const wchar_t* end);

#endif // BOOST_LOG_USE_WCHAR_T

} // namespace log

} // namespace boost

#endif // BOOST_LOG_NO_SETTINGS_PARSERS_SUPPORT
