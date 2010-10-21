/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/forces/lj.hpp>
#include <halmd/mdsim/gpu/forces/lj_kernel.hpp>
#include <halmd/mdsim/gpu/forces/pair_short_ranged_kernel.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

/**
 * Initialise Lennard-Jones potential parameters
 */
template <int dimension, typename float_type>
lj_potential<dimension, float_type>::lj_potential(
    unsigned ntype
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , r_cut_sigma_(ntype, ntype)
  , r_cut_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , sigma2_(ntype, ntype)
  , en_cut_(ntype, ntype)
  , g_param_(epsilon_.data().size())
{
    // FIXME support any number of types
    for (unsigned i = 0; i < std::min(ntype, 2U); ++i) {
        for (unsigned j = i; j < std::min(ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            r_cut_sigma_(i, j) = cutoff[i + j];
        }
    }

    // precalculate derived parameters
    for (unsigned i = 0; i < ntype; ++i) {
        for (unsigned j = i; j < ntype; ++j) {
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            // energy shift due to truncation at cutoff length
            float_type rri_cut = std::pow(r_cut_sigma_(i, j), -2);
            float_type r6i_cut = rri_cut * rri_cut * rri_cut;
            en_cut_(i, j) = 4 * epsilon_(i, j) * r6i_cut * (r6i_cut - 1);
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential core width: σ = " << sigma_);
    LOG("potential cutoff length: r_c = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);

    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[lj_kernel::EPSILON] = epsilon_.data()[i];
        p[lj_kernel::RR_CUT] = rr_cut_.data()[i];
        p[lj_kernel::SIGMA2] = sigma2_.data()[i];
        p[lj_kernel::EN_CUT] = en_cut_.data()[i];
        param[i] = p;
    }

    cuda::copy(param, g_param_);
}

template <int dimension, typename float_type>
lj<dimension, float_type>::lj(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
)
  // dependency injection
  : _Base(
        make_shared<potential_type>(particle->ntype, cutoff, epsilon, sigma)
      , particle
      , box
    )
{
}

template <int dimension, typename float_type>
void lj<dimension, float_type>::luaopen(lua_State* L)
{
    typedef typename _Base::_Base _Base2;
    typedef typename _Base2::_Base _Base3;
    using namespace luabind;
    string class_name("lj_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    namespace_("forces")
                    [
                        // skip auxiliary class _Base [_Base=forces::pair_short_ranged]
                        class_<lj, shared_ptr<_Base3>, bases<_Base2, _Base3> >(class_name.c_str())
                            .def(constructor<
                                shared_ptr<particle_type>
                              , shared_ptr<box_type>
                              , array<float, 3> const&
                              , array<float, 3> const&
                              , array<float, 3> const&
                            >())
                            .def("register_runtimes", &lj::register_runtimes)
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(2) //< distance of derived to base class
    [
        &lj<3, float>::luaopen
    ]
    [
        &lj<2, float>::luaopen
    ];
}

// explicit instantiation
template class lj<3, float>;
template class lj<2, float>;

}}} // namespace mdsim::gpu::forces

} // namespace halmd
