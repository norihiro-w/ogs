/**
 * \author Norihiro Watanabe
 * \date   2013-09-06
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <limits>
#include <valarray>
#include <algorithm>

#ifdef OGS_USE_EIGEN
#include <Eigen>
#endif

#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/FemEnums.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeData.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalCoordinatesMapping.h"

#include "../TestTools.h"

using namespace NumLib;

#ifdef OGS_USE_EIGEN
TEST(NumLib, FemNaturalCoordinatesMappingQuad4WithEigen)
{
    static const double eps = std::numeric_limits<double>::epsilon();

    // Eigen matrix types
    static const unsigned dim = 2;
    static const unsigned e_nnodes = 4;
    typedef Eigen::Matrix<double, e_nnodes, 1> NodalVector;
    typedef Eigen::Matrix<double, dim, e_nnodes, Eigen::RowMajor> DimNodalMatrix;
    typedef Eigen::Matrix<double, dim, dim, Eigen::RowMajor> DimMatrix;

    // Shape data type
    typedef ShapeData<NodalVector,DimNodalMatrix,DimMatrix> ShapeDataType;

    // Coordinates mapping type
    typedef NaturalCoordinatesMapping<MeshLib::Quad, ShapeQuad4, ShapeDataType> NaturalCoordinatesMappingType;

    // prepare a quad mesh element

    // identical to natural coordinates
    {
        MeshLib::Node** nodes = new MeshLib::Node*[4];
        MeshLib::Node node1( 1.0,  1.0, 0.0, 0);
        MeshLib::Node node2(-1.0,  1.0, 0.0, 1);
        MeshLib::Node node3(-1.0, -1.0, 0.0, 2);
        MeshLib::Node node4( 1.0, -1.0, 0.0, 3);
        nodes[0] = &node1;
        nodes[1] = &node2;
        nodes[2] = &node3;
        nodes[3] = &node4;
        MeshLib::Quad e(nodes);
        ShapeDataType shape(e.getDimension(), e.getNNodes());
        NaturalCoordinatesMappingType mapping(e);

        double r[2] = {0.5, 0.5};
        mapping.computeMappingMatrices(r, shape);
        double exp_N[]= {0.5625, 0.1875, 0.0625, 0.1875};
        double exp_dNdr[]= {0.375, -0.375, -0.125, 0.125, 0.375, 0.125, -0.125, -0.375};
        double exp_J[]= {1.0, 0.0, 0.0, 1.0};

        ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
        ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
        ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
        ASSERT_ARRAY_NEAR(exp_J, shape.invJ.data(), shape.invJ.size(), eps);
        ASSERT_NEAR(1.0, shape.detJ, eps);
        ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdx.data(), shape.dNdx.size(), eps);
    }

    // identical to natural coordinates (clockwise node ordering, which is invalid)
    {
        MeshLib::Node** nodes = new MeshLib::Node*[4];
        MeshLib::Node node1( 1.0,  1.0, 0.0, 0);
        MeshLib::Node node2(-1.0,  1.0, 0.0, 1);
        MeshLib::Node node3(-1.0, -1.0, 0.0, 2);
        MeshLib::Node node4( 1.0, -1.0, 0.0, 3);
        nodes[0] = &node1;
        nodes[1] = &node4;
        nodes[2] = &node3;
        nodes[3] = &node2;
        MeshLib::Quad e(nodes);
        ShapeDataType shape(e.getDimension(), e.getNNodes());
        NaturalCoordinatesMappingType mapping(e);

        double r[2] = {0.5, 0.5};
        mapping.computeMappingMatrices(r, shape);
        //std::cout << shape;
        double exp_N[]= {0.5625, 0.1875, 0.0625, 0.1875};
        double exp_dNdr[]= {0.375, -0.375, -0.125, 0.125, 0.375, 0.125, -0.125, -0.375};
        double exp_J[]= {0.0, 1.0, 1.0, 0.0};
        // Inverse of the Jacobian matrix doesn't exist
        double exp_invJ[dim*dim]= {0.0};
        double exp_dNdx[dim*e_nnodes]= {0.0};

        ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
        ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
        ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
        ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), eps);
        ASSERT_NEAR(-1.0, shape.detJ, eps);
        ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), eps);
    }

    // irregular shape
    {
        MeshLib::Node** nodes = new MeshLib::Node*[4];
        MeshLib::Node node1(-0.5, -0.5, 0.0, 0);
        MeshLib::Node node2( 0.6, -0.6, 0.0, 1);
        MeshLib::Node node3( 0.5,  0.4, 0.0, 2);
        MeshLib::Node node4(-0.3,  0.1, 0.0, 3);
        nodes[0] = &node1;
        nodes[1] = &node2;
        nodes[2] = &node3;
        nodes[3] = &node4;
        MeshLib::Quad e(nodes);
        ShapeDataType shape(e.getDimension(), e.getNNodes());
        NaturalCoordinatesMappingType mapping(e);

        double r[2] = {0.5, 0.5};
        mapping.computeMappingMatrices(r, shape);
//        std::cout << shape;
        double exp_N[]= {0.5625, 0.1875, 0.0625, 0.1875};
        double exp_dNdr[]= {0.375, -0.375, -0.125, 0.125, 0.375, 0.125, -0.125, -0.375};
        double exp_J[]= {-0.5125, 0.0, -0.0625, -0.35};
        double exp_invJ[]= {-1.9512195121951219, 0.0, 0.3484320557491290, -2.8571428571428572};
        double exp_dNdx[]= {-0.73170731707317072, 0.73170731707317072, 0.243902439024390, -0.24390243902439029, -0.940766550522648, -0.48780487804878048, 0.313588850174216, 1.1149825783972125};

        ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
        ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
        ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
        ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), eps);
        ASSERT_NEAR(0.179375, shape.detJ, eps);
        ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), eps);
    }



}
#endif

