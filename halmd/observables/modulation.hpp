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
  : public std::unary_function<fixed_vector<float_type, dimension>, float_type>
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
  : public std::unary_function<fixed_vector<float_type, dimension>, float_type>
{
public:
    exponential(float_type kappa, float_type offset, float_type box_height)
      : kappa_(kappa), offset_(offset), box_height_(box_height) {}

    HALMD_GPU_ENABLED float_type operator() (fixed_vector<float_type, dimension> const& r) const
    {
        float_type z = r[dimension-1];
#ifdef __CUDACC__
        z -= rintf(fdividef(z, box_height_)) * box_height_;     // fold into periodic box
        float_type arg = kappa_ * (z - offset_);
        return arg > 0 ? expf(-arg) : 1;
#else
        z -= rint(z / box_height_) * box_height_;
        float_type arg = kappa_ * (z - offset_);
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
    float_type box_height_;
};

/** exponential slab: @f$ f(\vec r) = [\exp(-\kappa (r_z - z_0)) - 1] \Theta(\kappa (r_z - z_0)) + 1
 *  \quad \mathrm{for} \quad \sgn(\kappa) (r_z - z_0) < w $@f
 */
template <int dimension, typename float_type>
class exponential_slab
  : public std::unary_function<fixed_vector<float_type, dimension>, float_type>
{
public:
    exponential_slab(float_type kappa, float_type width, float_type offset, float_type box_height)
      : kappa_(kappa), kappa_width_(std::abs(kappa * width)), offset_(offset), box_height_(box_height) {}

    HALMD_GPU_ENABLED float_type operator() (fixed_vector<float_type, dimension> const& r) const
    {
        float_type z = r[dimension-1];
#ifdef __CUDACC__
        z -= rintf(fdividef(z, box_height_)) * box_height_;     // fold into periodic box
        float_type arg = kappa_ * (z - offset_);
        return arg > 0 ? (arg < kappa_width_ ? expf(-arg) : 0) : 1;
#else
        z -= rint(z / box_height_) * box_height_;
        float_type arg = kappa_ * (z - offset_);
        return arg > 0 ? (arg < kappa_width_ ? exp(-arg) : 0) : 1;
#endif
    }

#ifndef __CUDACC__
    std::string message() const
    {
        std::ostringstream s;
        s << "use exponential slab modulation: κ = " << kappa_
          << ", w = " << kappa_width_ / kappa_ << ", z0 = " << offset_
          << ", precision = " << demangled_name<float_type>();
        return s.str();
    }
#endif

private:
    float_type kappa_;
    float_type kappa_width_;
    float_type offset_;
    float_type box_height_;
};

/** catenary: @f$ f(\vec r) = [\cosh(\kappa (r_z - z_0))) / \cosh(\kappa w / 2)) - 1] \Theta(w / 2 - |r_z - z_0|) + 1 $@f*/
template <int dimension, typename float_type>
class catenary
  : public std::unary_function<fixed_vector<float_type, dimension>, float_type>
{
public:
    catenary(float_type kappa, float_type width, float_type offset, float_type box_height)
      : kappa_(kappa), width_2_(width / 2), offset_(offset), box_height_(box_height)
      , norm_(cosh(kappa_ * width_2_))
    {}

    HALMD_GPU_ENABLED float_type operator() (fixed_vector<float_type, dimension> const& r) const
    {
        float_type z = r[dimension-1];
#ifdef __CUDACC__
        z -= rintf(fdividef(z, box_height_)) * box_height_;     // fold into periodic box
        z -= offset_;
        return fabsf(z) < width_2_ ? coshf(kappa_ * z) / norm_ : 1;
#else
        z -= rint(z / box_height_) * box_height_;
        z -= offset_;
        return std::abs(z) < width_2_ ? cosh(kappa_ * z) / norm_ : 1;
#endif
    }

#ifndef __CUDACC__
    std::string message() const
    {
        std::ostringstream s;
        s << "use catenary modulation: κ = " << kappa_ << ", w = " << 2 * width_2_ << ", z0 = " << offset_
          << ", precision = " << demangled_name<float_type>();
        return s.str();
    }
#endif

private:
    float_type kappa_;
    float_type width_2_;
    float_type offset_;
    float_type box_height_;
    float_type norm_;
};

} // namespace modulation
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_MODULATION_HPP */
