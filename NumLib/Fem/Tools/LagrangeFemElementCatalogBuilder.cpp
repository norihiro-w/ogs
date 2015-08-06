/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "LagrangeFemElementCatalogBuilder.h"

#include "NumLib/Fem/FiniteElement/FiniteElementType.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "FemElementCatalog.h"
#include "MeshElementShapeToFemElementType.h"

namespace NumLib
{

void LagrangeFemElementCatalogBuilder::construct(FemElementCatalog &feCatalog, MeshElementShapeToFemElementType &shape2fetype)
{
#ifdef OGS_USE_EIGEN
    feCatalog.registerFeType<NumLib::FeLINE2<EigenDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::LINE2);
    feCatalog.registerFeType<NumLib::FeQUAD4<EigenDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::QUAD4);
    feCatalog.registerFeType<NumLib::FeTRI3<EigenDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::TRI3);
    feCatalog.registerFeType<NumLib::FeTET4<EigenDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::TET4);
#else
    feCatalog.registerFeType<NumLib::FeLINE2<BlazeDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::LINE2);
    feCatalog.registerFeType<NumLib::FeQUAD4<BlazeDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::QUAD4);
    feCatalog.registerFeType<NumLib::FeTRI3<BlazeDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::TRI3);
    feCatalog.registerFeType<NumLib::FeTET4<BlazeDynamicShapeMatrixPolicy, 3>::type>(NumLib::FeType::TET4);
#endif
    shape2fetype.addFeType(MeshLib::MeshElemType::LINE, 1, FeType::LINE2);
    shape2fetype.addFeType(MeshLib::MeshElemType::QUAD, 1, FeType::QUAD4);
    shape2fetype.addFeType(MeshLib::MeshElemType::TRIANGLE, 1, FeType::TRI3);
    shape2fetype.addFeType(MeshLib::MeshElemType::TETRAHEDRON, 1, FeType::TET4);
}

}
