/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef SHAPEDATA_H_
#define SHAPEDATA_H_

#include "../FemEnums.h"

namespace NumLib
{

/**
 * \brief Coordinate mapping matrices at particular location
 *
 * Mapping is done between physical coordinates (x,y,z) and natural coordinates (r,s,t)
 *
 * \tparam T_N      Vector type for shape functions
 * \tparam T_DN     Matrix type for gradient of shape functions
 * \tparam T_J      Jacobian matrix type
 */
template <class T_N, class T_DN, class T_J>
struct ShapeData
{
    typedef T_N ShapeType;
    typedef T_DN DShapeType;
    typedef T_J JacobianType;

    /// shape function N(r)
    ShapeType N;

    /// gradient of shape functions, dN(r)/dr
    DShapeType dNdr;

    /// gradient of shape functions, dN(r)/dx
    DShapeType dNdx;

    /// Jacobian matrix, J=dx/dr
    JacobianType J;

    /// inverse of the Jacobian
    JacobianType invJ;

    /// determinant of the Jacobian
    double detJ;

    /**
     * Initialize matrices and vectors
     *
     * @param dim       Spatial dimension
     * @param n_nodes   The number of element nodes
     */
    ShapeData(std::size_t dim, std::size_t n_nodes)
    : N(n_nodes), dNdr(dim, n_nodes), dNdx(dim, n_nodes), J(dim, dim),
      invJ(dim, dim), detJ(.0) {}

    ~ShapeData() {}

    /**
     * reset data with zero
     *
     * @param fields    bit flags to specify fields to be initialized
     */
    void setZero(const ShapeFieldType fields = SHAPE_ALL);

    /**
     * writes the matrix entries into the output stream
     * @param out the output stream
     */
    void write (std::ostream& out) const;

private:
    template<class T>
    void setZero(T &mat)
    {
        mat.setZero(mat.rows(), mat.cols());
    }

    void setZero(ShapeType &vec)
    {
        vec.setZero(vec.size());
    }

};

template <class T_N, class T_DN, class T_J>
void ShapeData<T_N, T_DN, T_J>::setZero(const ShapeFieldType fields)
{
    if (fields & SHAPE_N)
        setZero(N);
    if (!(fields & SHAPE_DNDX)) return;
    setZero(dNdr);
    setZero(dNdx);
    setZero(J);
    setZero(invJ);
    detJ = .0;
}

template <class T_N, class T_DN, class T_J>
void ShapeData<T_N, T_DN, T_J>::write(std::ostream& out) const
{
    out << "N   :\n" << N << "\n";
    out << "dNdr:\n" << dNdr << "\n";
    out << "J   :\n" << J << "\n";
    out << "invJ:\n" << invJ << "\n";
    out << "|J| : " << detJ << "\n";
    out << "dNdx:\n" << dNdx << "\n";
}

template <class T_N, class T_DN, class T_J>
std::ostream& operator<< (std::ostream &os, const ShapeData<T_N, T_DN, T_J> &shape)
{
    shape.write (os);
    return os;
}

} // NumLib

#endif //SHAPEDATA_H_
