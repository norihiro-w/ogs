//============================================================================
// Name        : fracture_map.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

struct RasterData
{
    std::size_t ncols;
    std::size_t nrows;
    double xllcorner;
    double yllcorner;
    double cellsize;
    double NODATA_Value;
    std::vector<double> values;
};



void outputQuadMesh(const RasterData &raster, const std::string &filename)
{
    const std::size_t n_org_cells = raster.ncols * raster.nrows;
    const std::size_t n_new_nodes = (raster.ncols+1)*(raster.nrows+1);
    const std::size_t n_new_eles = n_org_cells;

    std::ofstream os(filename.c_str());
    os.precision(std::numeric_limits<double>::digits10);
    os << "#FEM_MSH\n";
    os << " $PCS_TYPE\n";
    os << "  NO_PCS\n";
    os << " $NODES\n";
    os << "  " << n_new_nodes << "\n";
    const double dL = raster.cellsize;
    const double offset_x = raster.xllcorner;
    const double offset_y = raster.yllcorner+raster.nrows*dL;
    // from top left to bottom right
    std::size_t node_id = 0;
    for (std::size_t i=0; i<raster.nrows+1; i++) {
        for (std::size_t j=0; j<raster.ncols+1; j++) {
            os << "  " << node_id++ << "\t" << j*dL+offset_x << "\t" <<  offset_y - ((double)i)*dL << "\t" << 0 << "\n";
        }
    }
    os << " $ELEMENTS\n";
    os << "  " << n_new_eles << "\n";
    std::size_t ele_id = 0;
    for (std::size_t i=0; i<raster.nrows; i++) {
        for (std::size_t j=0; j<raster.ncols; j++) {
            os << "  " << ele_id++ << " 0 quad "  << i*(raster.ncols+1) + j << "\t" <<  (i+1)*(raster.ncols+1) + j << "\t" << (i+1)*(raster.ncols+1) + (j+1) << "\t" << i*(raster.ncols+1) + (j+1) << "\n";
        }
    }
    os << "#STOP\n";
    os.close();
}

void outputElementProperty(const RasterData &raster, const std::string &filename)
{
    std::ofstream os(filename.c_str());
    os.precision(std::numeric_limits<double>::digits10);
    os << "#MEDIUM_PROPERTIES_DISTRIBUTED\n";
    os << " $MSH_TYPE\n  NO_PCS\n";
    os << " $MMP_TYPE\n  APERTURE\n";
    os << " $DIS_TYPE\n  ELEMENT\n";
    os << " $DATA\n";
    for (size_t i=0; i<raster.values.size(); i++)
        os << "  " << i << " " << raster.values[i] << "\n";
    os << "#STOP\n";
    os.close();
}

void outputVTK(const RasterData &raster, const std::string &filename)
{
    const std::size_t n_org_cells = raster.ncols * raster.nrows;
    const std::size_t n_new_nodes = (raster.ncols+1)*(raster.nrows+1);

    std::ofstream vtk_file(filename.c_str());
    vtk_file.precision(std::numeric_limits<double>::digits10);
    vtk_file << "# vtk DataFile Version 2.0\n";
    vtk_file << "Unstructured Grid\n";
    vtk_file << "ASCII\n\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";
    vtk_file << "POINTS "<< n_new_nodes << " float\n";
    const double dL = raster.cellsize;
    const double offset_x = raster.xllcorner;
    const double offset_y = raster.yllcorner+raster.nrows*dL;
    // from top left to bottom right
    for (std::size_t i=0; i<raster.nrows+1; i++) {
        for (std::size_t j=0; j<raster.ncols+1; j++) {
            vtk_file << "  " << j*dL+offset_x << " " <<  offset_y - ((double)i)*dL << " " << 0 << "\n";
        }
    }
    vtk_file << "\n";

    // write element header
    vtk_file << "CELLS " << n_org_cells << " " << n_org_cells*(1+4) << "\n";

    // write elements
    for (std::size_t i=0; i<raster.nrows; i++) {
        for (std::size_t j=0; j<raster.ncols; j++) {
            vtk_file << "4  " << i*(raster.ncols+1) + j << " " <<  (i+1)*(raster.ncols+1) + j << " " << (i+1)*(raster.ncols+1) + (j+1) << " " << i*(raster.ncols+1) + (j+1) << "\n";
        }
    }
    vtk_file << "\n";
    vtk_file << "CELL_TYPES " <<n_org_cells << "\n";
    for (size_t i=0; i<raster.nrows*raster.ncols; i++) {
       vtk_file  << "9 \n";
    }
    vtk_file << "\n";

    vtk_file << "CELL_DATA " <<n_org_cells << "\n";
    vtk_file << "SCALARS " << "ELE_ATTRIBUTE" << " float 1" << "\n";
    vtk_file << "LOOKUP_TABLE default" <<"\n";
    //....................................................................
    for (size_t i=0; i<raster.nrows*raster.ncols; i++) {
        vtk_file <<" "<< raster.values[i] << "\n";
    }

    vtk_file.close();
}

void readRaster(const std::string &filename, RasterData &raster)
{
    std::ifstream is(filename.c_str());
    std::string dummy;
    is >> dummy >> raster.ncols;
    is >> dummy >> raster.nrows;
    is >> dummy >> raster.xllcorner;
    is >> dummy >> raster.yllcorner;
    is >> dummy >> raster.cellsize;
    is >> dummy >> raster.NODATA_Value;
    const std::size_t n_data = raster.ncols * raster.nrows;
    raster.values.resize(n_data);
    for (std::size_t i=0; i<n_data; i++)
        is >> raster.values[i];

    is.close();
}

void convertMMtoM(RasterData &data)
{
    const double fac = 1e-3;
    data.xllcorner *= fac;
    data.yllcorner *= fac;
    data.cellsize *= fac;
    for (size_t i=0; i<data.values.size(); i++)
        data.values[i] *= fac;
}

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cout << "USAGE:" << std::endl;
        std::cout << "\t" << argv[0] << " (input ASCII file) (output mesh file) (output prop file) (output vtk file)" << std::endl;
        return 0;
    }
    std::string in_filename(argv[1]);
    std::string out_filename(argv[2]);
    std::string out_prop_filename(argv[3]);
    std::string out_vtk_filename(argv[4]);

    RasterData data;
    std::cout << "-> reading a raster file..." << std::endl;
    readRaster(in_filename, data);
    std::cout << "-> converting a unit from mm to m..." << std::endl;
    convertMMtoM(data);
    std::cout << "-> creating a quad mesh file..." << std::endl;
    outputQuadMesh(data, out_filename);
    std::cout << "-> creating a element property file..." << std::endl;
    outputElementProperty(data, out_prop_filename);
    std::cout << "-> creating a VTK file..." << std::endl;
    outputVTK(data, out_vtk_filename);

	return 0;
}
