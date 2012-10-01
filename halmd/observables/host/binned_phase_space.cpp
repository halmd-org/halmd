/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/bind.hpp>

#include <halmd/observables/host/binned_phase_space.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
binned_phase_space<dimension, float_type>::binned_phase_space(
    shared_ptr<sample_type> sample
  , shared_ptr<box_type const> box
  , shared_ptr<clock_type const> clock
  , shared_ptr<logger_type> logger
)
  // module dependencies
  : sample_(sample)
  , clock_(clock)
  , logger_(logger)
{
    vector_type box_length = static_cast<vector_type>(box->length());
    sample_->cell_length = element_div(box_length, static_cast<vector_type>(sample_->nbin));
    sample_->cell_origin = static_cast<vector_type>(box->origin());

    // set up grid of cell positions axis-wise
    for (unsigned axis = 0; axis < dimension; ++axis) {
        double spacing = sample_->cell_length[axis];
        double offset = sample_->cell_origin[axis];
        position_[axis].reserve(sample_->nbin[axis]);
        for (unsigned i = 0; i < sample_->nbin[axis]; ++i) {
            position_[axis].push_back((i + 0.5) * spacing + offset);
        }
    }
}

/**
 * Bin phase space sample into spatial cells
 */
template <int dimension, typename float_type>
void binned_phase_space<dimension, float_type>::acquire()
{
    if (sample_->step == clock_->step()) {
        LOG_TRACE("sample is up to date");
        return;
    }

    // trigger update of phase space sample
    on_acquire_();

    LOG_TRACE("acquire sample");
    scoped_timer_type timer(runtime_.acquire);

    // empty cell lists without memory reallocation
    cell_array_type& cell_array = sample_->cell_array;
    for_each(
        cell_array.data(), cell_array.data() + cell_array.num_elements()
      , bind(&cell_type::clear, _1)
    );

    typedef typename data_sample_type::sample_vector sample_vector;
    sample_vector const& r_sample = sample_->position_data();

    // add particles to cells
    vector_type const& cell_origin = sample_->cell_origin;
    vector_type const& cell_length = sample_->cell_length;
    cell_size_type const& nbin = sample_->nbin;
    for (size_t i = 0; i < r_sample.size(); ++i) {
        vector_type const& r = r_sample[i];
        // compute cell index from position
        //
        // Since the positions are extended beyond the periodic simulation box,
        // we have to account for modulo operations on large negative values;
        // in C99, "a % n" yields values in (-n, n)
        cell_diff_type index = element_mod(
            static_cast<cell_diff_type>(floor(element_div(r - cell_origin, cell_length)))
          , static_cast<cell_diff_type>(nbin)
        );
        for (int j = 0; j < dimension; ++j) {
            if (index[j] < 0) { index[j] += nbin[j]; }
        }
        cell_array(static_cast<cell_size_type>(index)).push_back(i);
    }

    sample_->step = clock_->step();
}

template <int dimension, typename float_type>
void binned_phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("binned_phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                class_<binned_phase_space, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("acquire", &runtime::acquire)
                    ]
                    .def_readonly("runtime", &binned_phase_space::runtime_)
            ]

          , def("binned_phase_space", &make_shared<binned_phase_space
                , shared_ptr<sample_type>
                , shared_ptr<box_type const>
                , shared_ptr<clock_type const>
                , shared_ptr<logger_type>
            >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_binned_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    binned_phase_space<3, double>::luaopen(L);
    binned_phase_space<2, double>::luaopen(L);
#endif
    binned_phase_space<3, float>::luaopen(L);
    binned_phase_space<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class binned_phase_space<3, double>;
template class binned_phase_space<2, double>;
#endif
template class binned_phase_space<3, float>;
template class binned_phase_space<2, float>;

} // namespace host
} // namespace observables
} // namespace halmd
