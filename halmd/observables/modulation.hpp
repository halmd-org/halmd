/*
 * Copyright © 2013 Felix Höfling
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

#ifndef HALMD_OBSERVABLES_MODULATION_HPP
#define HALMD_OBSERVABLES_MODULATION_HPP

#include <functional>
#ifndef __CUDACC__
#   include <sstream>
#endif

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#ifndef __CUDACC__
#   include <halmd/utility/demangle.hpp>
#endif

namespace halmd {
namespace observables {
namespace modulation {

/**
 * define modulation factors for computation of density modes
 */

// typedef BOOST_TYPEOF((boost::lambda::_1, float_type(1))) unity;

/** unity: @f$ f(\vec r) = 1 $@f*/
template <int dimension, typename float_type>
struct unity
  : public std::unary_function<float_type, fixed_vector<float_type, dimension> const&>
{
    HALMD_GPU_ENABLED float_type operator() (fixed_vector<float_type, dimension> const& r) const
    {
        return 1;
    }

#ifndef __CUDACC__
    std::string message() const
    {
        return std::string();
    }
#endif
};

/** exponential: @f$ f(\vec r) = [\exp(-\kappa (r_z - z_0)) - 1] \Theta(\kappa (r_z - z_0)) + 1 $@f*/
template <int dimension, typename float_type>
class exponential
  : public std::unary_function<float_type, fixed_vector<float_type, dimension> const&>
{
public:
    exponential(float_type kappa, float_type offset)
      : kappa_(kappa), offset_(offset) {}

    HALMD_GPU_ENABLED float_type operator() (fixed_vector<float_type, dimension> const& r) const
    {
        float_type arg = kappa_ * (r[dimension-1] - offset_);
#ifdef __CUDACC__
        return arg > 0 ? expf(-arg) : 1;
#else
        return arg > 0 ? exp(-arg) : 1;
#endif
    }

#ifndef __CUDACC__
    std::string message() const
    {
        std::ostringstream s;
        s << "use exponential modulation: κ = " << kappa_ << ", z0 = " << offset_
          << ", precision = " << demangled_name<float_type>();
        return s.str();
    }
#endif

private:
    float_type kappa_;
    float_type offset_;
};

/** catenary: @f$ f(\vec r) = [\cosh(\kappa (r_z - z_0))) / \cosh(\kappa (w - z_0)) - 1] \Theta(w - |r_z - z_0|) + 1 $@f*/
template <int dimension, typename float_type>
class catenary
  : public std::unary_function<float_type, fixed_vector<float_type, dimension> const&>
{
public:
    catenary(float_type kappa, float_type width, float_type offset)
      : kappa_(kappa), width_(width), offset_(offset), norm_(cosh(kappa_ * width_)) {}

    HALMD_GPU_ENABLED float_type operator() (fixed_vector<float_type, dimension> const& r) const
    {
        float_type z = r[dimension-1] - offset_;
#ifdef __CUDACC__
        return fabsf(z) < width_ ? coshf(kappa_ * z) / norm_ : 1;
#else
        return std::abs(z) < width_ ? cosh(kappa_ * z) / norm_ : 1;
#endif
    }

#ifndef __CUDACC__
    std::string message() const
    {
        std::ostringstream s;
        s << "use catenary modulation: κ = " << kappa_ << ", w = " << width_ << ", z0 = " << offset_
          << ", precision = " << demangled_name<float_type>();
        return s.str();
    }
#endif

private:
    float_type kappa_;
    float_type width_;
    float_type offset_;
    float_type norm_;
};

} // namespace modulation
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_MODULATION_HPP */
