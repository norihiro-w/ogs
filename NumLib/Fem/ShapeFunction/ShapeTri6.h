/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPETRI6_H_
#define SHAPETRI6_H_

#include "MeshLib/Elements/Tri.h"

namespace NumLib
{

/**
 *  Shape function for a 6-nodes triangle element
 *
 */
class ShapeTri6
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] N    a vector of calculated shape function.
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X &r, T_N &N);

    /**
     * Evaluate derivatives of the shape function at the given point
     * The point coordinates in \c r are not used.
     *
     * @param [in]  r    point coordinates
     * @param [out] dN   a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X &r, T_N &dN);

    using MeshElement = MeshLib::Tri6;
    static const std::size_t DIM = MeshElement::dimension;
    static const std::size_t NPOINTS = MeshElement::n_all_nodes;
};

}

#include "ShapeTri6-impl.h"

#endif //SHAPETRI6_H_
