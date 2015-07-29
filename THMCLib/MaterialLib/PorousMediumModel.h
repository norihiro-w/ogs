/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PorousMedia.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Elements/Element.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "IMedium.h"
#include "BasicModel.h"
#include "MaterialProperties.h"
#include "PorousMediumThermalConductivity.h"
#include "PorousMediumHeatCapacity.h"

namespace NumLib
{
class ITXFunction;
}

namespace MaterialLib
{

typedef BasicModelTensor HydraulicConductivity;
typedef BasicModelTensor Permeability;
typedef TemplateBasicModel<double> Porosity;
typedef TemplateBasicModel<double> SpecificStorage;
typedef TemplateBasicModel<double> GeometricArea;
typedef TemplateBasicModel<double> MassDispersivityLong;
typedef TemplateBasicModel<double> MassDispersivityTrans;


struct PorousMediumModel : public IMedium
{
    HydraulicConductivity* hydraulic_conductivity;
    Permeability* permeability;
    Porosity* porosity;
    SpecificStorage* storage;
    GeometricArea* geo_area;
    MassDispersivityLong* dispersivity_long;
    MassDispersivityTrans* dispersivity_trans;
    PorousMediumHeatCapacity* heat_capacity;
    PorousMediumThermalConductivity* thermal_conductivity;

    PorousMediumModel() : hydraulic_conductivity(nullptr), permeability(nullptr), porosity(nullptr),
            storage(nullptr), geo_area(nullptr), dispersivity_long(nullptr), dispersivity_trans(nullptr),
            heat_capacity(nullptr), thermal_conductivity(nullptr)
    {
    }

    virtual ~PorousMediumModel()
    {
        delete hydraulic_conductivity;
        delete permeability;
        delete porosity;
        delete storage;
        delete geo_area;
        delete dispersivity_long;
        delete dispersivity_trans;
        delete heat_capacity;
        delete thermal_conductivity;
    }

    virtual MediumType getMediumType() const {return MediumType::PorousMedium;};

    PorousMediumProperty operator()(const StateVariables &var, unsigned global_dim, const MeshLib::Element &e, const SolidProperty &s, const FluidProperty &f) const
    {
        std::vector<FluidProperty> vec_f(1, f);
        return (*this)(var, global_dim, e, s, vec_f);
    }

    PorousMediumProperty operator()(const StateVariables &var, unsigned global_dim, const MeshLib::Element &e, const SolidProperty &s, const std::vector<FluidProperty> &f) const
    {
        const MathLib::RotationMatrix* matR = (global_dim==e.getDimension()) ? nullptr : &e.getMappedLocalCoordinates().getRotationMatrixToGlobal();
        PorousMediumProperty pm;
        if (geo_area) pm.geo_area = (*geo_area)(&var);
        if (permeability) pm.k = (*permeability)(&var, e.getDimension(), global_dim, matR);
        if (porosity) pm.n = (*porosity)(&var);
        if (storage) pm.Ss = (*storage)(&var);
        if (heat_capacity) pm.Cp = (*heat_capacity)(&var, pm, s, f);
        if (thermal_conductivity) pm.lamba = (*thermal_conductivity)(&var, pm, s, f, global_dim);
        return pm;
    }
};

} //end
 
