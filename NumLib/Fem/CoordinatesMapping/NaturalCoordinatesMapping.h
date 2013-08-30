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


#ifndef NATURALCOORDINATESMAPPING_H_
#define NATURALCOORDINATESMAPPING_H_


#include "ShapeData.h"
#include "../FemEnums.h"

namespace NumLib
{

/**
 * \brief Mapping element shapes to natural coordinates
 *
 * This class is responsible for mapping element shapes in physical coordinates
 * (x,y,z) to those in natural coordinates (r,s,t).
 * - Given physical coordinates should correspond to dimensions of the element, 
 *   i.e. (x,y) for triangles
 * - _nodes_coords(r,s,t) = N(r,s,t) x_i
 *
 * @tparam T_MESH_ELEMENT   Mesh element class
 * @tparam T_SHAPE_FUNC     Shape function class
 * @tparam T_SHAPE_DATA     Shape data class
 */
template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
class NaturalCoordinatesMapping
{
public:
    typedef T_MESH_ELEMENT MeshElementType;
    typedef T_SHAPE_DATA ShapeDataType;

    NaturalCoordinatesMapping() : _ele(nullptr) {}

    /// initialize this object with the given mesh element
    explicit NaturalCoordinatesMapping(const MeshElementType &ele);

    ~NaturalCoordinatesMapping() {}

    /// reset a mesh element
    void reset(const MeshElementType &ele);

    /// compute mapping matrices at the given location in natural coordinates
    ///
    /// @param natural_pt   Location in natural coordinates
    /// @param shape        Shape data where calculated shape functions are stored
    /// @param fields       Bit flags to specify fields to be calculated
    void computeMappingMatrices(const double* natural_pt, ShapeDataType &shape, const ShapeFieldType fields=SHAPE_ALL) const;

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{_nodes_coords} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param shape        Shape data evaluated at the natural coordinates
    /// @param physical_pt  Location in physical coordinates
    void mapToPhysicalCoordinates(const ShapeDataType &shape, double* physical_pt) const;

    /// compute natural coordinates at the given physical coordinates.
    /// Assuming \f$ r=0 \f$ at \f$ _nodes_coords = \bar{_nodes_coords}_{avg} \f$, natural coordinates can be calculated as
    /// \f[
    ///    \mathbf{r} = (\mathbf{J}^-1)^T * (\mathbf{_nodes_coords} - \bar{\mathbf{_nodes_coords}}_{avg})
    /// \f]
    ///
    /// @param shape        Shape data evaluated at some location
    /// @param physical_pt  Location in physical coordinates
    /// @param natural_pt   Location in natural coordinates
    void mapToNaturalCoordinates(const ShapeDataType &shape, const double* physical_pt, double* natural_pt) const;

private:
    const MeshElementType* _ele;
};

} //end namespace

#include "NaturalCoordinatesMapping.tpp"


#endif //NATURALCOORDINATESMAPPING_H_
