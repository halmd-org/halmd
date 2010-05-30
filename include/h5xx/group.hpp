/*
 * Copyright © 2010  Peter Colberg
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

#ifndef H5XX_GROUP_HPP
#define H5XX_GROUP_HPP

#include <boost/iterator/iterator_facade.hpp>
#include <string>
#include <vector>

#include <h5xx/id.hpp>

namespace h5xx
{

/**
 * HDF5 group base class for file or group
 */
class basic_group
  : public virtual basic_id
{
public:
    explicit basic_group(hid_t hid)
      : basic_id(hid)
    {}
};

class basic_group_iterator;

/**
 * This class wraps a HDF5 group.
 */
class group
  : public virtual basic_id
  , public basic_group
{
public:
    typedef basic_group_iterator iterator;

    /**
     * open or create group
     */
    group(basic_group const& parent, std::string const& path)
      : basic_id(open(parent, path))
      , basic_group(basic_id::hid())
    {}

    /**
     * open or create group within group
     */
    group& operator/=(std::string const& path)
    {
        group child(*this, path);
        basic_id::swap(child);
        return *this;
    }

    /**
     * check whether a group exists
     */
    static bool exists(basic_group const& parent, std::string const& path)
    {
        hid_t hid;
        try {
#ifdef H5XX_USE_16_API
            H5XX_CHECK(hid = H5Gopen(parent.hid(), path.c_str()));
#else
            H5XX_CHECK(hid = H5Gopen(parent.hid(), path.c_str(), H5P_DEFAULT));
#endif
        }
        catch (error const& e)
        {
            if (!e.count(std::make_pair(H5E_SYM, H5E_NOTFOUND))) {
                throw;
            }
            return false;
        }
        basic_group child(hid); // wrap identifier for destruction
        return true;
    }

private:
    /**
     * open or create group
     */
    static basic_id open(basic_id const& parent, std::string const& path)
    {
        hid_t hid;
        try {
#ifdef H5XX_USE_16_API
            H5XX_CHECK(hid = H5Gopen(parent.hid(), path.c_str()));
#else
            H5XX_CHECK(hid = H5Gopen(parent.hid(), path.c_str(), H5P_DEFAULT));
#endif
        }
        catch (error const& e)
        {
            if (!e.count(std::make_pair(H5E_SYM, H5E_NOTFOUND))) {
                throw;
            }
#ifdef H5XX_USE_16_API
            H5XX_CHECK(hid = H5Gcreate(parent.hid(), path.c_str(), 0));
#else
            H5XX_CHECK(hid = H5Gcreate(parent.hid(), path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
#endif
        }
        return basic_id(hid);
    }
};

/**
 * open or create group within file
 */
inline group operator/(basic_group const& parent, std::string const& path)
{
    return group(parent, path);
}

/**
 * HDF5 group iterator
 */
class basic_group_iterator
  : public boost::iterator_facade<
        basic_group_iterator
      , group
      , boost::forward_traversal_tag
      , group
    >
{
public:
    /**
     * start iterator
     */
    basic_group_iterator(basic_group const& parent_)
      : parent_(&parent_)
      , pos_(seek(0)) // first group object in sequence
    {}

    /**
     * end iterator
     */
    basic_group_iterator()
      : parent_(NULL)
      , pos_(-1) // undefined
    {}

private:
    friend class boost::iterator_core_access;

    /**
     * increment iterator and seek to next group object
     */
    void increment()
    {
        pos_ = seek(++pos_);
    }

    /**
     * seek to next group object starting at given index
     */
    hsize_t seek(hsize_t pos) const
    {
        while (!end(pos)) {
#ifdef H5XX_USE_16_API
            int ret;
            H5XX_CHECK(ret = H5Gget_objtype_by_idx(parent_->hid(), pos));
            if (H5G_GROUP == ret) {
                break;
            }
#else
            H5O_info_t info;
            H5XX_CHECK(H5Oget_info_by_idx(parent_->hid(), ".", H5_INDEX_NAME, H5_ITER_INC, pos, &info, H5P_DEFAULT));
            if (H5O_TYPE_GROUP == info.type) {
                break;
            }
#endif
            ++pos;
        }
        return pos;
    }

    /*
     * check if position is end of sequence
     */
    bool end(hsize_t pos) const
    {
        if (parent_ != NULL) {
#ifdef H5XX_USE_16_API
            hsize_t count;
            H5XX_CHECK(H5Gget_num_objs(parent_->hid(), &count));
            return pos == count;
#else
            H5G_info_t info;
            H5XX_CHECK(H5Gget_info(parent_->hid(), &info));
            return pos == info.nlinks;
#endif
        }
        return true;
    }

    /**
     * Check equality with another iterator.
     *
     * This is used to compare against the end iterator, which is
     * equal to the default constructed iterator.
     */
    bool equal(basic_group_iterator const& other) const
    {
        if (end(pos_)) {
            return other.end(other.pos_);
        }
        if (other.end(other.pos_)) {
            return false;
        }
        return *parent_ == *other.parent_ && pos_ == other.pos_;
    }

    /**
     * This does not actually dereference a pointer to a group
     * instance, but instead constructs a group instance given
     * the parent group and the name of the child group at the
     * iterator's sequence position.
     *
     * Be sure to dereference the iterator with (*it).…, and not
     * it->…, as the second form causes the invokation of the
     * basic_id copy constructor and is therefore less efficient.
     */
    group dereference() const
    {
        ssize_t size; // excludes NULL terminator
#ifdef H5XX_USE_16_API
        H5XX_CHECK(size = H5Gget_objname_by_idx(parent_->hid(), pos_, NULL, 0));
#else
        H5XX_CHECK(size = H5Lget_name_by_idx(parent_->hid(), ".", H5_INDEX_NAME, H5_ITER_INC, pos_, NULL, 0, H5P_DEFAULT));
#endif
        std::vector<char> name(size + 1); // includes NULL terminator
#ifdef H5XX_USE_16_API
        H5XX_CHECK(H5Gget_objname_by_idx(parent_->hid(), pos_, name.data(), name.size()));
#else
        H5XX_CHECK(H5Lget_name_by_idx(parent_->hid(), ".", H5_INDEX_NAME, H5_ITER_INC, pos_, name.data(), name.size(), H5P_DEFAULT));
#endif
        return *parent_ / name.data();
    }

    /** parent group object to iterate over */
    basic_group const* parent_;
    /** position within alphabetically sorted index sequence */
    hsize_t pos_;
};

} // namespace h5xx

#endif /* ! H5XX_GROUP_HPP */
