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


#include <iostream>
#include <cassert>

#include "logog/include/logog.hpp"

#include "MathLib/LinAlg/VectorNorms.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

namespace NumLib
{

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::NaturalCoordinatesMapping(const MeshElementType &ele)
: _ele(&ele)
{
    this->reset(ele);
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::reset(const MeshElementType &ele)
{
    _ele = &ele;
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::computeMappingMatrices(const double* natural_pt, ShapeDataType &prop, const ShapeFieldType fields) const
{
    //prepare
    const std::size_t dim = _ele->getDimension();
    const std::size_t nnodes = _ele->getNNodes();
    prop.setZero(fields);

    //shape, dshape/dr
    if ( fields & SHAPE_N )
        T_SHAPE_FUNC::computeShapeFunction(natural_pt, prop.N);

    if ( !(fields & SHAPE_DNDX)) return;

    double* dNdr = prop.dNdr.data();
    T_SHAPE_FUNC::computeGradShapeFunction(natural_pt, dNdr);

    //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
    for (std::size_t k=0; k<nnodes; k++) {
        const double* xyz = _ele->getNode(k)->getCoords();
        // outer product of dNdr and xyz for a particular node
        for (std::size_t j_x=0; j_x<dim; j_x++) {
            for (std::size_t i_r=0; i_r<dim; i_r++) {
                prop.J(i_r,j_x) += prop.dNdr(i_r,k) * xyz[j_x];
            }
        }
    }

    //|J|, J^-1, dshape/dx
    prop.detJ = prop.J.determinant();
    if (prop.detJ>.0) {
        prop.invJ = prop.J.inverse();
        prop.dNdx = prop.invJ * prop.dNdr;
    }
#ifndef NDEBUG
    if (prop.detJ<=.0)
        ERR("***error: det|J|=%e is not positive.\n", prop.detJ);
#endif
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::mapToPhysicalCoordinates(const ShapeDataType &prop, double* physical_pt) const
{
    const std::size_t dim = _ele->getDimension();
    const std::size_t nnodes = _ele->getNNodes();

    for (std::size_t i=0; i<dim; i++) {
        physical_pt[i] = .0;
        for (std::size_t j=0; j<nnodes; j++)
            physical_pt[i] += prop.N(j) * _ele->getNode(j)->getCoords()[i];
    }
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::mapToNaturalCoordinates(const ShapeDataType &prop, const double* physical_pt, double* natural_pt) const
{
    const std::size_t dim = _ele->getDimension();
    const std::size_t nnodes = _ele->getNNodes();

#ifdef OGS5_APPROACH
    // calculate dx which is relative coordinates from the element center
    std::vector<double> center_pt(dim, .0);
    // x_avg = sum_i {x_i} / n
    for (std::size_t i=0; i<dim; i++)
        for (std::size_t j=0; j<nnodes; j++)
            center_pt[i] += _ele->getNode(j)->getCoords()[i];
    for (std::size_t i=0; i<dim; i++)
        center_pt[i] /= (double)nnodes;
    for (auto v : center_pt) std::cout << v << " "; std::cout << std::endl;
    // dx = pt - x_avg
    std::vector<double> dx(dim, .0);
    for (std::size_t i=0; i<dim; i++)
        dx[i] = physical_pt[i] - center_pt[i];
    for (auto v : dx) std::cout << v << " "; std::cout << std::endl;

    // r = invJ^T * dx
    for (std::size_t i=0; i<dim; i++)
        natural_pt[i] = 0.0;
    for (std::size_t i=0; i<dim; i++) {
        for (std::size_t j=0; j<dim; j++)
            natural_pt[i] += prop.invJ(j, i) * dx[j];
    }

#else
    // solve R(r)=N(r)*{X} - x = 0 with Newton-Raphson
    std::cout << "*** Newton iteration start" << std::endl;
    std::valarray<double> N(nnodes);
    MathLib::DenseMatrix<double> dN(dim, nnodes), J(dim, dim);
    std::valarray<double> r(dim), dr(dim), R(dim);
    const std::size_t max_iterations = 10;
    const double eps = std::numeric_limits<double>::epsilon();
    for (std::size_t i=0; i<max_iterations; i++) {
        // evaluate residual
        R = .0;
        T_SHAPE_FUNC::computeShapeFunction(r, N);
        for (std::size_t j=0; j<dim; j++) {
            for (std::size_t k=0; k<nnodes; k++)
                R[j] += N[k]*_ele->getNode(k)->getCoords()[j];
            R[j] -= physical_pt[j];
        }
        double error = MathLib::normEuklid(&R[0], 2);
        std::cout << i << ": error = " << error << std::endl;
        if (error<eps)
            break;

        // Jacobain = dNdr(r)*{X}
        double* pdN = &dN(0,0);
        T_SHAPE_FUNC::computeGradShapeFunction(r, pdN);
        for (std::size_t j=0; j<dim; j++) {
            for (std::size_t k=0; k<dim; k++) {
                J(j,k) = .0;
                for (std::size_t l=0; l<nnodes; l++)
                    J(j,k) += dN(k,l)*_ele->getNode(l)->getCoords()[j];
            }
        }
        MathLib::GaussAlgorithm<MathLib::DenseMatrix<double>, std::valarray<double>> gauss(J);
        gauss.solve(R, dr);
        r -= dr;
        std::cout << "r= "; for (auto v: r) std::cout << v << " "; std::cout << std::endl;
    }
    for (std::size_t i=0; i<dim; i++)
        natural_pt[i] = r[i];
#endif
}

} //end
