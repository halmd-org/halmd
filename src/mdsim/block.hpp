/* Block algorithm parameters
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_BLOCK_HPP
#define MDSIM_BLOCK_HPP

#include <algorithm>
#include <cmath>
#include "exception.hpp"
#include "log.hpp"


namespace mdsim {

/**
 * Block algorithm parameters
 */
template <unsigned dimension, typename T>
class block_param
{
public:
    /** set total number of simulation steps */
    void steps(uint64_t const& value, float const& timestep);
    /** set total simulation time */
    void time(float const& value, float const& timestep);
    /** set block size */
    void block_size(unsigned int const& value);
    /** set maximum number of samples per block */
    void max_samples(uint64_t const& value);

    /** returns total number of simulation steps */
    uint64_t const& steps() const;
    /** returns total simulation time */
    float const& time() const;
    /** returns block size */
    unsigned int const& block_size() const;
    /** returns block shift */
    unsigned int const& block_shift() const;
    /** returns block count */
    unsigned int const& block_count() const;
    /** returns maximum number of samples per block */
    uint64_t const& max_samples() const;

    /** returns simulation time belonging to block sample at given block offset */
    float timegrid(unsigned int block, unsigned int offset) const;

private:
    /** simulation timestep */
    float timestep_;
    /** total number of simulation steps */
    uint64_t steps_;
    /** total simulation time */
    float time_;
    /** block size */
    unsigned int block_size_;
    /** block shift */
    unsigned int block_shift_;
    /** block count */
    unsigned int block_count_;
    /** maximum number of samples per block */
    uint64_t max_samples_;
};


/**
 * set total number of simulation steps
 */
template <unsigned dimension, typename T>
void block_param<dimension, T>::steps(uint64_t const& value, float const& timestep)
{
    // set total number of simulation steps
    steps_ = value;
    LOG("total number of simulation steps: " << steps_);
    // set simulation timestep
    timestep_ = timestep;
    // derive total simulation time
    time_ = value * timestep_;
    LOG("total simulation time: " << time_);
}

/**
 * set total simulation time
 */
template <unsigned dimension, typename T>
void block_param<dimension, T>::time(float const& value, float const& timestep)
{
    // set total simulation time
    time_ = value;
    LOG("total simulation time: " << time_);
    // set simulation timestep
    timestep_ = timestep;
    // derive total number of simulation steps
    steps_ = std::ceil(time_ / timestep_);
    LOG("total number of simulation steps: " << steps_);
}

/**
 * set block size
 */
template <unsigned dimension, typename T>
void block_param<dimension, T>::block_size(unsigned int const& value)
{
    // set block size
    block_size_ = value;
    LOG("block size: " << block_size_);

    // derive block shift from block size
    block_shift_ = std::floor(std::sqrt(block_size_));
    LOG("block shift: " << block_shift_);
    if (block_shift_ < 2) {
	throw exception("computed block shift is less than 2, larger block size required");
    }

    // derive block count from block size and block shift
    block_count_ = 0;
    for (unsigned int n = block_size_; n <= steps_; n *= block_size_) {
	block_count_++;
	if ((n * block_shift_) > steps_) {
	    break;
	}
	block_count_ ++;
    }
    LOG("block count: " << block_count_);
    if (block_count_ < 2) {
	throw exception("computed block count is less than 2, more simulations steps required");
    }
}

/**
 * set maximum number of samples per block
 */
template <unsigned dimension, typename T>
void block_param<dimension, T>::max_samples(uint64_t const& value)
{
    max_samples_ = std::min(value, steps_);
    if (max_samples_ != value) {
	LOG_WARNING("overriding maximum number of samples with number of simulation steps");
    }
    LOG("maximum number of samples per block: " << max_samples_);
    if (max_samples_ < block_size_) {
	throw exception("maximum number of samples must not be smaller than block size");
    }
}

/**
 * returns total number of simulation steps
 */
template <unsigned dimension, typename T>
uint64_t const& block_param<dimension, T>::steps() const
{
    return steps_;
}

/**
 * returns total simulation time
 */
template <unsigned dimension, typename T>
float const& block_param<dimension, T>::time() const
{
    return time_;
}

/**
 * returns block size
 */
template <unsigned dimension, typename T>
unsigned int const& block_param<dimension, T>::block_size() const
{
    return block_size_;
}

/**
 * returns block shift
 */
template <unsigned dimension, typename T>
unsigned int const& block_param<dimension, T>::block_shift() const
{
    return block_shift_;
}

/**
 * returns block count
 */
template <unsigned dimension, typename T>
unsigned int const& block_param<dimension, T>::block_count() const
{
    return block_count_;
}

/**
 * returns maximum number of samples per block
 */
template <unsigned dimension, typename T>
uint64_t const& block_param<dimension, T>::max_samples() const
{
    return max_samples_;
}

/**
 * returns simulation time belonging to block sample at given block offset
 */
template <unsigned dimension, typename T>
float block_param<dimension, T>::timegrid(unsigned int block, unsigned int offset) const
{
    return timestep_ * std::pow(block_size_, block / 2.f) * (offset + 1) * ((block % 2) ? block_shift_ : 1) ;
}

} // namespace mdsim

#endif /* ! MDSIM_BLOCK_HPP */
