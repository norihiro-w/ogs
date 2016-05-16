/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_COMPONENTGLOBALINDEXDICT_H_
#define ASSEMBLERLIB_COMPONENTGLOBALINDEXDICT_H_

#include <iostream>
#include <limits>

#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

#include "MeshLib/Location.h"

namespace AssemblerLib
{

/// \internal
namespace detail
{

template <typename IndexType>
struct Line
{
    MeshLib::Location location;

    // Physical component
    std::size_t comp_id;

    // Position in global matrix or vector
    IndexType global_index;

    Line(MeshLib::Location const& l, std::size_t c, IndexType i)
    : location(l), comp_id(c), global_index(i)
    {}

    Line(MeshLib::Location const& l, std::size_t c)
    : location(l), comp_id(c),
        global_index(std::numeric_limits<IndexType>::max())
    {}

    explicit Line(MeshLib::Location const& l)
    : location(l),
        comp_id(std::numeric_limits<std::size_t>::max()),
        global_index(std::numeric_limits<IndexType>::max())
    {}

    friend std::ostream& operator<<(std::ostream& os, Line const& l)
    {
        return os << l.location << ", " << l.comp_id << ", " << l.global_index;
    }
};

template <typename IndexType>
struct LineByLocationComparator
{
    bool operator()(Line<IndexType> const& a, Line<IndexType> const& b) const
    {
        return a.location < b.location;
    }
};

template <typename IndexType>
struct LineByLocationAndComponentComparator
{
    bool operator()(Line<IndexType> const& a, Line<IndexType> const& b) const
    {
        if (a.location < b.location)
            return true;
        if (b.location < a.location)
            return false;

        // a.loc == b.loc
        return a.comp_id < b.comp_id;
    }
};

struct ByLocation {};
struct ByLocationAndComponent {};
struct ByComponent {};
struct ByGlobalIndex {};

template <typename IndexType>
using ComponentGlobalIndexDict = boost::multi_index::multi_index_container<
        Line<IndexType>,
        boost::multi_index::indexed_by
        <
            boost::multi_index::ordered_unique
            <
                boost::multi_index::tag<ByLocationAndComponent>,
                boost::multi_index::identity<Line<IndexType>>,
                LineByLocationAndComponentComparator<IndexType>
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByLocation>,
                boost::multi_index::identity<Line<IndexType>>,
                LineByLocationComparator<IndexType>
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByComponent>,
                boost::multi_index::member<Line<IndexType>, std::size_t, &Line<IndexType>::comp_id>
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByGlobalIndex>,
                boost::multi_index::member<Line<IndexType>, IndexType, &Line<IndexType>::global_index>
            >
        >
    >;

}    // namespace detail
}    // namespace AssemblerLib

#endif  // ASSEMBLERLIB_COMPONENTGLOBALINDEXDICT_H_
