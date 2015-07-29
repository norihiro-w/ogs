# Example modules for THMC in porous media

## Process simulation modules
- GROUNDWATER_FLOW
  - hydraulic head based formulation of groundawter flow in confined aquifer
  - Input: -
  - Output: Hydraulic head
- LIQUID_FLOW
  - liquid phase pressure based formulation of groundwater flow in confined aquifer
  - Input: -
  - Output: Liquid phase pressure
- HEAT_TRANSPORT
  - Heat transport in saturated porous media
  - Input: Mean fluid velocity
  - Output: Temperature
- MASS_TRANSPORT
  - Mass transport in saturated porous media
  - Input: Mean fluid velocity
  - Output: Concentration
- DEFORMATION
  - Displacement based formulation of linear elastic deformation (total deformation)
  - Input: -
  - Output: Displacement
- INCREMENTAL_DEFORMATION
  - Incremental displacement based formulation of linear elastic deformation (incremental deformation)
  - Input: -
  - Output: Displacement, Strain, Stress
- DEFORMATION_FLOW
  - Monolithic formulation of poroelasticity (fully saturated case, total deformation, Darcy flow)
  - Input: -
  - Output: Displacement, Pressure

## Post processing modules
- HEAD_TO_ELEMENT_VELOCITY
  - Calculation of Darcy velocity at integration points from hydraulic head
  - Input: Hydraulic head
  - Output: Darcy velocity at integration points
- PRESSURE_TO_ELEMENT_VELOCITY
  - Calculation of Darcy velocity at integration points from liquid-phase pressure
  - Input: Liquid phase pressure
  - Output: Darcy velocity at integration points
- PRESSURE_TO_HEAD
  - Calculation of hydraulic head from liquid-phase pressure
  - Input: Liquid phase pressure
  - Output: Hydraulic head
- ELEMENT_STRESS_STRAIN
  - Calculation of strain and stress at integration points from displacement
  - Input: Displacement
  - Output: Stress and strain at integration points
- NODAL_STRESS_STRAIN
  - Extrapolation of strain and stress at integration points to those at nodes
  - Input: Stress and strain at integration points
  - Output: Stress and strain at nodes