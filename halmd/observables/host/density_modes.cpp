/*
 * Copyright © 2011  Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#include <boost/foreach.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/observables/host/density_modes.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using boost::fusion::at_key;
using namespace std;

namespace halmd
{
namespace observables { namespace host
{

template <int dimension, typename float_type>
density_modes<dimension, float_type>::density_modes(
    shared_ptr<density_modes_sample_type> rho_sample
  , shared_ptr<trajectory_sample_type> trajectory_sample
  , shared_ptr<trajectory_stream_type> trajectory_stream
  , vector<double> const& wavenumbers
  , vector_type const& box_length
  , double tolerance
  , unsigned int max_count
)
    // dependency injection
  : rho_sample(rho_sample)
  , trajectory_sample(trajectory_sample)
  , trajectory_stream(trajectory_stream)
    // member initialisation
  , wavevectors_(wavenumbers, box_length, tolerance, max_count)
{
    // number of particle types must agree
    assert(rho_sample->rho.size() == trajectory_sample->r.size());

    // allocate memory
    unsigned int ntype = rho_sample->rho.size();
    unsigned int nq = wavevectors_.values().size();
    for (unsigned int i = 0; i < ntype; ++i) {
        typedef typename density_modes_sample_type::mode_vector_type mode_vector_type;
        rho_sample->rho[i].reset(new mode_vector_type(nq));
    }
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void density_modes<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * register data request from a sink
 */
template <int dimension, typename float_type>
void density_modes<dimension, float_type>::register_request(uint64_t step, function<void(uint64_t)> callback)
{
    LOG_TRACE("[density_modes] request received for step " << step);
    request_.insert(make_pair(step, callback));

    // in order to satisfy the request, we need to issue requests to our sources,
    // but only if not already done for this timestamp 'step'
    if (issued_request_.find(step) == issued_request_.end()) {
        LOG_TRACE("[density_modes] issue requests for step " << step);
        trajectory_stream->register_request(step, bind(&density_modes::notify, this, _1));
        issued_request_.insert(step);
    }
}

/**
 * notification when trajectory data are available
 *
 */
template <int dimension, typename float_type>
void density_modes<dimension, float_type>::notify(uint64_t step)
{
    LOG_TRACE("[density_modes] notification in step " << step);

    // remove from list of open requests
    assert(issued_request_.find(step) != issued_request_.end());
    issued_request_.erase(step);

    // transform trajectory data to density modes
    acquire(step);  // FIXME rename to transform() or compute()

    // find requests with matching timestamp 'step'
    typedef request_container_type::iterator iterator_type;
    std::pair<iterator_type, iterator_type> range = request_.equal_range(step);

    // notify all callbacks of range
    for (iterator_type it = range.first; it != range.second; ++it) {
        it->second(step);
    }
    // remove fulfilled requests from list
    request_.erase(range.first, range.second);
}

/**
 * Acquire sample of all density modes from trajectory sample
 */
template <int dimension, typename float_type>
void density_modes<dimension, float_type>::acquire(double time)
{
    scoped_timer<timer> timer_(at_key<sample_>(runtime_));

    // do nothing if we're up to date
    if (rho_sample->time == time) return;

    typedef typename trajectory_sample_type::sample_vector_ptr positions_vector_ptr_type;
    typedef typename density_modes_sample_type::mode_vector_type mode_vector_type;

    // trigger update of trajectory sample
    if (trajectory_sample->time != time) {
        // FIXME
    }

    // compute density modes separately for each particle type
    // 1st loop: iterate over particle types
    unsigned int type = 0;
    BOOST_FOREACH (positions_vector_ptr_type const r_sample, trajectory_sample->r) {
        mode_vector_type& rho_vector = *rho_sample->rho[type]; //< dereference shared_ptr
        // initialise result array
        fill(rho_vector.begin(), rho_vector.end(), 0);
        // compute sum of exponentials: rho_q = sum_r exp(-i q·r)
        // 2nd loop: iterate over particles of the same type
        BOOST_FOREACH (vector_type const& r, *r_sample) {
            typename mode_vector_type::iterator rho_q = rho_vector.begin();
            typedef pair<double, vector_type> map_value_type; // pair: (wavenumber, wavevector)
            // 3rd loop: iterate over wavevectors
            BOOST_FOREACH (map_value_type const& q_pair, wavevectors_.values()) {
                float_type q_r = inner_prod(static_cast<vector_type>(q_pair.second), r);
                *rho_q++ += mode_type(cos(q_r), -sin(q_r));
            }
        }
        ++type;
    }
    rho_sample->time = time;
}

template <int dimension, typename float_type>
void density_modes<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("density_modes_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                namespace_("host")
                [
                    class_<density_modes, shared_ptr<_Base>, _Base>(class_name.c_str())
                        .def(constructor<
                            shared_ptr<density_modes_sample_type>
                          , shared_ptr<trajectory_sample_type>
                          , shared_ptr<trajectory_stream_type>
                          , vector<double> const&
                          , vector_type const&, double, unsigned int
                        >())
                        .def("register_runtimes", &density_modes::register_runtimes)
                        .property("tolerance", &density_modes::tolerance)
                        .property("maximum_count", &density_modes::maximum_count)
                ]
            ]
        ]
    ];
}

namespace  // limit symbols to translation unit
{

__attribute__ ((constructor)) void register_lua()
{
    lua_wrapper::register_(1)	//< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &density_modes<3, double>::luaopen
    ]
    [
        &density_modes<2, double>::luaopen
    ];
#else
    [
        &density_modes<3, float>::luaopen
    ]
    [
        &density_modes<2, float>::luaopen
    ];
#endif
}

}  // namespace

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class density_modes<3, double>;
template class density_modes<2, double>;
#else
template class density_modes<3, float>;
template class density_modes<2, float>;
#endif

}}  // namespace observables::host

}  // namespace halmd
