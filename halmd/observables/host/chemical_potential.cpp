/*
 * Copyright © 2016 Felix Höfling
 * Copyright © 2016 Arthur Straube
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cmath>
#include <sstream>

#include <halmd/observables/host/chemical_potential.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace host {

template <int dimension, typename float_type>
chemical_potential<dimension, float_type>::chemical_potential(
    std::shared_ptr<particle_type> test_particle
  , std::shared_ptr<position_type> position
  , double temperature
  , std::vector<unsigned int> ntest_particle
  , std::shared_ptr<halmd::logger> logger
)
  : test_particle_(test_particle)
  , position_(position)
  , ntest_particle_(ntest_particle)
  , mu_ex_(ntest_particle.size())
  , logger_(logger)
{
    set_temperature(temperature);

    std::ostringstream str;
    std::copy(ntest_particle_.begin(), ntest_particle_.end(), std::ostream_iterator<double>(str, " "));
    LOG("number of test particles: " << str.str());

    // construct particle instance for test particles
    unsigned int npart = std::accumulate(ntest_particle_.begin(), ntest_particle_.end(), 0);
    if (npart != test_particle_->nparticle()) {
        throw std::invalid_argument("test particle instance has mismatching number of particles");
    }

    // assign particle positions
    set_position();

    // assign particle species' according to values in ntest_particle_
    auto species = make_cache_mutable(test_particle_->species());
    auto it = species->begin();
    unsigned int s = 0;
    for (unsigned int n : ntest_particle_) {     // iterate over species counts
        std::fill(it, it + n, s);
        it += n;
        ++s;
    }
    assert(it == species->end());
}

template <int dimension, typename float_type>
void chemical_potential<dimension, float_type>::set_temperature(double temperature)
{
    temperature_ = temperature;
    LOG("temperature for computation of chemical potential: " << temperature_);
}

// FIXME support particle groups
template <int dimension, typename float_type>
typename chemical_potential<dimension, float_type>::result_type const&
chemical_potential<dimension, float_type>::sample()
{
    // check cache of potential_energy, this may trigger a recalculation of the energy
    cache<en_pot_array_type> const& en_pot_cache = test_particle_->potential_energy();

    if (mu_ex_cache_ != en_pot_cache) {
        LOG_TRACE("sample excess chemical potential");
        scoped_timer_type timer(runtime_.sample);

        // iterate over energies of test particles, for each species separately
        auto const& en_pot = read_cache(en_pot_cache);
        auto it = en_pot.begin();

        for (unsigned int s = 0; s < ntest_particle_.size(); ++s) {
            unsigned int n = ntest_particle_[s];           // number of particles of species s

            // accumulate exp(-E/kT) for particles in range (it, it + n)
            accumulator<double> acc;
            std::for_each(it, it + n, [&](double x) {
                if (std::isfinite(x)) {         // catch NaN/Inf etc.
                    acc(exp(-x / temperature_));
                }
                else {
                    acc(0);
                }
            });
            it += n;

            // compute μ_ex = -kT log( <exp(-E/kT)> ) and store with statistics
            double Z = mean(acc);
            mu_ex_[s] = accumulator<double>(
                - temperature_ * log(Z)                                     // mean
              , pow(temperature_ / Z, 2) * variance(acc)                    // variance
              , count(acc)                                                  // number of test particles of this species
            );
        }

        mu_ex_cache_ = en_pot_cache;
    }
    return mu_ex_;
}

template <typename chemical_potential_type>
static std::function<accumulator<double> const& ()>
wrap_sample(std::shared_ptr<chemical_potential_type> self, unsigned int species)
{
    if (species >= self->result_size()) {
        LOG_ERROR("requested species too large: " << species);
        throw std::invalid_argument("index exceeds size of result array");
    }
    return [=]() -> accumulator<double> const& {
        return self->sample()[species];
    };
}

template <int dimension, typename float_type>
void chemical_potential<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<chemical_potential>()
                .def("sample", &wrap_sample<chemical_potential>)
                .def("set_position", &chemical_potential::set_position)
                .property("temperature", &chemical_potential::temperature, &chemical_potential::set_temperature)
                .property("test_particle", &chemical_potential::test_particle)
                .property("result_size", &chemical_potential::result_size)
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("sample", &runtime::sample)
                ]
                .def_readonly("runtime", &chemical_potential::runtime_)

          , def("chemical_potential", &std::make_shared<chemical_potential
                , std::shared_ptr<particle_type>
                , std::shared_ptr<position_type>
                , double
                , std::vector<unsigned int>
                , std::shared_ptr<logger>
             >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_chemical_potential(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    chemical_potential<3, double>::luaopen(L);
    chemical_potential<2, double>::luaopen(L);
#else
    chemical_potential<3, float>::luaopen(L);
    chemical_potential<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class chemical_potential<3, double>;
template class chemical_potential<2, double>;
#else
template class chemical_potential<3, float>;
template class chemical_potential<2, float>;
#endif

} // namespace host
} // namespace observables
} // namespace halmd
