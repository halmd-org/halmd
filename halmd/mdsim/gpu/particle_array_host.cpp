/*
 * Copyright © 2017 Daniel Kirchner
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

#include <halmd/mdsim/gpu/particle_array_host.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template<typename T>
particle_array_host<T>::particle_array_host(
  std::shared_ptr<particle_array_gpu_base> const& parent
, size_t offset
, size_t stride
, bool coalesced) : parent_(parent), offset_(offset), stride_(stride), coalesced_(coalesced) {
}

template<typename T>
particle_array_host<T>::~particle_array_host() {
}

template<typename T>
std::shared_ptr<particle_array_host<T>> particle_array_host<T>::cast(std::shared_ptr<particle_array_host_base> base) {
    if(base->type() != typeid(T)) {
        throw std::runtime_error("invalid cast");
    }
    return std::static_pointer_cast<particle_array_host>(base);
}

template<typename T>
void particle_array_host<T>::set_lua(luaponte::object table)
{
    auto mem = coalesced_ ? parent_->get_host_data() : parent_->get_host_memory();
    size_t offset = offset_;
    // TODO: handle size mismatch & ghost particles
    for(luaponte::iterator it(table), end; it != end && offset < mem.size(); ++it) {
        helper::set(mem, offset, luaponte::object_cast<T>(*it));
        offset += stride_;
    }
    parent_->set_host_data(mem);
}

template<typename T>
luaponte::object particle_array_host<T>::get_lua(lua_State *L) const
{
    auto data = parent_->get_host_data();
    luaponte::object table = luaponte::newtable(L);
    std::size_t offset = offset_;
    // TODO: handle ghost particles
    for(std::size_t j = 1; j <= parent_->nparticle(); j++) {
        auto&& value = helper::get(data, offset);
        table[j++] = boost::cref(value);
        offset += stride_;
    }
    return table;
}


template class particle_array_host<float>;
template class particle_array_host<fixed_vector<float, 2>>;
template class particle_array_host<fixed_vector<float, 3>>;
template class particle_array_host<fixed_vector<float, 4>>;
template class particle_array_host<unsigned int>;
template class particle_array_host<fixed_vector<unsigned int, 2>>;
template class particle_array_host<fixed_vector<unsigned int, 3>>;
template class particle_array_host<fixed_vector<unsigned int, 4>>;
template class particle_array_host<int>;
template class particle_array_host<fixed_vector<int, 2>>;
template class particle_array_host<fixed_vector<int, 3>>;
template class particle_array_host<fixed_vector< int, 4>>;
template class particle_array_host<stress_tensor_wrapper<typename type_traits<2, float>::stress_tensor_type>>;
template class particle_array_host<stress_tensor_wrapper<typename type_traits<3, float>::stress_tensor_type>>;

} // namespace gpu
} // namespace mdsim
} // namespace halmd
