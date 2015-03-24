/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Implementation of the MshLayerMapper class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "LogogSimpleFormatter.h"

// GeoLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"
#include "MeshSurfaceExtraction.h"
#include "MathTools.h"
#include "MeshEnums.h"

#include "readMeshFromFile.h"
#include "Legacy/MeshIO.h"
#include "FileIO/VtkIO/VtuInterface.h"

void addHalfLayers(	const std::vector<MeshLib::Node*> &surface_nodes, const std::vector<MeshLib::Element*> &surface_elems,
					unsigned nHalfLayers, double z_offset, double dz, double maxZ, unsigned &mat_id, bool mat_id_each_layer,
					std::vector<MeshLib::Node*> &new_nodes, std::vector<MeshLib::Element*> &new_elems, size_t &node_offset, size_t &elem_offset)
{
	const size_t nSurfaceNodes = surface_nodes.size();
	const size_t nSurfaceElems = surface_elems.size();

	for (unsigned layer_id = 0; layer_id <= nHalfLayers; ++layer_id)
	{
		if (layer_id==0 && node_offset > 0) continue;

		// add surface_nodes for new layer
		if (layer_id>0) z_offset += dz;
		for (unsigned i = 0; i < nSurfaceNodes; ++i)
		{
			const double* coords = surface_nodes[i]->getCoords();
			double new_z = coords[2]-z_offset;
			if (layer_id>0 && layer_id==nHalfLayers) {
				new_z = maxZ;
				if (dz<0) new_z += -dz;
			}
			new_nodes[node_offset+i] = new MeshLib::Node(coords[0], coords[1], new_z, node_offset+i);
		}

		// starting with 2nd layer create prism or hex elements connecting the last layer with the current one
		if (layer_id > 0)
		{
			size_t temp_node_offset = node_offset - nSurfaceNodes;
//			const unsigned mat_id (0);
//				const unsigned mat_id (nHalfLayers - layer_id);
//			std::cout << "-> add layer " << layer_id << " with mat_id " << mat_id << "\n";

			for (unsigned i = 0; i < nSurfaceElems; ++i)
			{
				const MeshLib::Element* sfc_elem( surface_elems[i] );
				if (sfc_elem->getDimension() == 2)
				{
					const unsigned nElemNodes(sfc_elem->getNNodes());
					MeshLib::Node** e_nodes = new MeshLib::Node*[2*nElemNodes];

					for (unsigned j=0; j<nElemNodes; ++j)
					{
						const unsigned node_id = sfc_elem->getNode(j)->getID() + temp_node_offset;
						if (dz > .0) {
							e_nodes[j] = new_nodes[node_id+nSurfaceNodes];
							e_nodes[j+nElemNodes] = new_nodes[node_id];
						} else {
							e_nodes[j+nElemNodes] = new_nodes[node_id+nSurfaceNodes];
							e_nodes[j] = new_nodes[node_id];
						}
					}
					if (sfc_elem->getGeomType() == MeshElemType::TRIANGLE)	// extrude triangles to prism
						new_elems[elem_offset+i] = new MeshLib::Prism(e_nodes, mat_id);
					else if (sfc_elem->getGeomType() == MeshElemType::QUAD)	// extrude quads to hexes
						new_elems[elem_offset+i] = new MeshLib::Hex(e_nodes, mat_id);
				}
				else
				{
					std::cout << "Warning in MshLayerMapper::CreateLayers() - Method can only handle 2D mesh elements ..." << std::endl;
					std::cout << "Skipping Element " << i << " of type \"" << MeshElemType2String(sfc_elem->getGeomType()) << "\"." << std::endl;
				}
			}
			elem_offset += nSurfaceElems;
			if (mat_id_each_layer) mat_id++;
		}
		node_offset += nSurfaceNodes;
//		std::cout << "-> node_offset = " << node_offset << ", elem_offset = " << elem_offset << "\n";
	}
}

void extrudeLayer( const std::vector<MeshLib::Node*> &surface_nodes, const std::vector<MeshLib::Element*> &surface_elems,
                    unsigned nHalfLayers, double z_offset, double dz, unsigned &mat_id, bool mat_id_each_layer,
                    std::vector<MeshLib::Node*> &new_nodes, std::vector<MeshLib::Element*> &new_elems, size_t &node_offset, size_t &elem_offset)
{
    const size_t nSurfaceNodes = surface_nodes.size();
    const size_t nSurfaceElems = surface_elems.size();

    for (unsigned layer_id = 0; layer_id <= nHalfLayers; ++layer_id)
    {
        if (layer_id==0 && node_offset > 0) continue;

        // add surface_nodes for new layer
        if (layer_id>0) z_offset += dz;
        for (unsigned i = 0; i < nSurfaceNodes; ++i)
        {
            const double* coords = surface_nodes[i]->getCoords();
            double new_z = coords[2]-z_offset;
            new_nodes[node_offset+i] = new MeshLib::Node(coords[0], coords[1], new_z, node_offset+i);
        }

        // starting with 2nd layer create prism or hex elements connecting the last layer with the current one
        if (layer_id > 0)
        {
            size_t temp_node_offset = node_offset - nSurfaceNodes;
//          const unsigned mat_id (0);
//              const unsigned mat_id (nHalfLayers - layer_id);
//          std::cout << "-> add layer " << layer_id << " with mat_id " << mat_id << "\n";

            for (unsigned i = 0; i < nSurfaceElems; ++i)
            {
                const MeshLib::Element* sfc_elem( surface_elems[i] );
                if (sfc_elem->getDimension() == 2)
                {
                    const unsigned nElemNodes(sfc_elem->getNNodes());
                    MeshLib::Node** e_nodes = new MeshLib::Node*[2*nElemNodes];

                    for (unsigned j=0; j<nElemNodes; ++j)
                    {
                        const unsigned node_id = sfc_elem->getNode(j)->getID() + temp_node_offset;
                        if (dz > .0) { // to down
                            e_nodes[j] = new_nodes[node_id+nSurfaceNodes];
                            e_nodes[j+nElemNodes] = new_nodes[node_id];
                        } else { // to up
                            e_nodes[j+nElemNodes] = new_nodes[node_id+nSurfaceNodes];
                            e_nodes[j] = new_nodes[node_id];
                        }
                    }
                    if (sfc_elem->getGeomType() == MeshElemType::TRIANGLE) {   // extrude triangles to prism
                        new_elems[elem_offset+i] = new MeshLib::Prism(e_nodes, mat_id);
                    } else if (sfc_elem->getGeomType() == MeshElemType::QUAD) { // extrude quads to hexes
                        new_elems[elem_offset+i] = new MeshLib::Hex(e_nodes, mat_id);
                    }
                }
                else
                {
                    std::cout << "Warning in MshLayerMapper::CreateLayers() - Method can only handle 2D mesh elements ..." << std::endl;
                    std::cout << "Skipping Element " << i << " of type \"" << MeshElemType2String(sfc_elem->getGeomType()) << "\"." << std::endl;
                }
            }
            elem_offset += nSurfaceElems;
            if (mat_id_each_layer) mat_id++;
        }
        node_offset += nSurfaceNodes;
//      std::cout << "-> node_offset = " << node_offset << ", elem_offset = " << elem_offset << "\n";
    }
}

void connect(	const std::vector<MeshLib::Node*> &lower_nodes, const std::vector<MeshLib::Element*> &lower_elems,
				const std::vector<MeshLib::Node*> &upper_nodes, const std::vector<MeshLib::Element*> &upper_elems,
					double dz, unsigned &mat_id,
					std::vector<MeshLib::Node*> &new_nodes, std::vector<MeshLib::Element*> &new_elems,
					size_t &node_offset, size_t &elem_offset)
{
//	std::cout << "-> add an intermediate layer\n";

	const size_t nUpperNodes = upper_nodes.size();
	const size_t nUpperElems = upper_elems.size();

	for (unsigned i = 0; i < nUpperNodes; ++i)
	{
		const double* coords = upper_nodes[i]->getCoords();
		new_nodes[node_offset+i] = new MeshLib::Node(coords[0], coords[1], coords[2] + dz, node_offset+i);
	}

	for (unsigned i = 0; i < nUpperElems; ++i)
	{
		const MeshLib::Element* sfc_elem( upper_elems[i] );
		if (sfc_elem->getDimension() == 2)
		{
			const unsigned nElemNodes(sfc_elem->getNNodes());
			MeshLib::Node** e_nodes = new MeshLib::Node*[2*nElemNodes];

			for (unsigned j=0; j<nElemNodes; ++j)
			{
				const unsigned node_id = sfc_elem->getNode(j)->getID() + node_offset;
				e_nodes[j+nElemNodes] = new_nodes[node_id];
				e_nodes[j] = new_nodes[sfc_elem->getNode(j)->getID()]; // lower surface nodes
			}
			if (sfc_elem->getGeomType() == MeshElemType::TRIANGLE)	// extrude triangles to prism
				new_elems[elem_offset+i] = new MeshLib::Prism(e_nodes, mat_id);
			else if (sfc_elem->getGeomType() == MeshElemType::QUAD)	// extrude quads to hexes
				new_elems[elem_offset+i] = new MeshLib::Hex(e_nodes, mat_id);
		}
		else
		{
			std::cout << "Warning in MshLayerMapper::CreateLayers() - Method can only handle 2D mesh elements ..." << std::endl;
			std::cout << "Skipping Element " << i << " of type \"" << MeshElemType2String(sfc_elem->getGeomType()) << "\"." << std::endl;
		}
	}

	node_offset += nUpperNodes;
	elem_offset += nUpperElems;
	mat_id++;
//	std::cout << "-> node_offset = " << node_offset << ", elem_offset = " << elem_offset << "\n";
}

void createVoidLayer(
				const std::vector<MeshLib::Node*> &lower_nodes, const std::vector<MeshLib::Element*> &lower_elems,
				const std::vector<MeshLib::Node*> &upper_nodes, const std::vector<MeshLib::Element*> &upper_elems,
					double dz, unsigned &mat_id,
					std::vector<MeshLib::Node*> &new_nodes, std::vector<MeshLib::Element*> &new_elems,
					size_t &node_offset, size_t &elem_offset)
{
//	std::cout << "-> add an intermediate layer\n";

	const size_t nSurfaceNodes = upper_nodes.size();
	const size_t nSurfaceElems = upper_elems.size();

	// add nodes, first lower nodes and then upper nodes
	for (unsigned i = 0; i < nSurfaceNodes; ++i)
	{
		const double* lower_coords = lower_nodes[i]->getCoords();
		new_nodes[node_offset+i] = new MeshLib::Node(lower_coords[0], lower_coords[1], lower_coords[2], node_offset+i);
	}
	node_offset += nSurfaceNodes;
	// add upper nodes + offset
	for (unsigned i = 0; i < nSurfaceNodes; ++i)
	{
		const double* upper_coords = upper_nodes[i]->getCoords();
		new_nodes[node_offset+i] = new MeshLib::Node(upper_coords[0], upper_coords[1], upper_coords[2] + dz, node_offset+i);
	}
	node_offset += nSurfaceNodes;

	// create elements
	for (unsigned i = 0; i < nSurfaceElems; ++i)
	{
		const MeshLib::Element* sfc_elem( upper_elems[i] );
		if (sfc_elem->getDimension() == 2)
		{
			const unsigned nElemNodes(sfc_elem->getNNodes());
			MeshLib::Node** e_nodes = new MeshLib::Node*[2*nElemNodes];
			// collect element nodes
			for (unsigned j=0; j<nElemNodes; ++j)
			{
				const unsigned node_id = sfc_elem->getNode(j)->getID() + nSurfaceNodes;
				// lower nodes first
				e_nodes[j] = new_nodes[sfc_elem->getNode(j)->getID()];
				// upper nodes second
				e_nodes[j+nElemNodes] = new_nodes[node_id];
			}
			// extrude surface elements to volume
			if (sfc_elem->getGeomType() == MeshElemType::TRIANGLE)
				new_elems[elem_offset+i] = new MeshLib::Prism(e_nodes, mat_id);
			else if (sfc_elem->getGeomType() == MeshElemType::QUAD)
				new_elems[elem_offset+i] = new MeshLib::Hex(e_nodes, mat_id);
		}
		else
		{
			std::cout << "Warning in MshLayerMapper::CreateLayers() - Method can only handle 2D mesh elements ..." << std::endl;
			std::cout << "Skipping Element " << i << " of type \"" << MeshElemType2String(sfc_elem->getGeomType()) << "\"." << std::endl;
		}
	}

	elem_offset += nSurfaceElems;
	mat_id++;
//	std::cout << "-> node_offset = " << node_offset << ", elem_offset = " << elem_offset << "\n";
}

MeshLib::Mesh* CreateLayers(const MeshLib::Mesh* uppesr_mesh, const MeshLib::Mesh* lower_mesh, double extZ, double dz)
{
	std::size_t nHalfLayers = extZ / dz;

	const size_t nUpperNodes = uppesr_mesh->getNNodes();
	const size_t nUpperElems = uppesr_mesh->getNElements();
	const std::vector<MeshLib::Node*> upper_nodes = uppesr_mesh->getNodes();
	const std::vector<MeshLib::Element*> upper_elems = uppesr_mesh->getElements();
	const size_t nLowerNodes = lower_mesh->getNNodes();
	const size_t nLowerElems = lower_mesh->getNElements();
	const std::vector<MeshLib::Node*> lower_nodes = lower_mesh->getNodes();
	const std::vector<MeshLib::Element*> lower_elems = lower_mesh->getElements();

	if (nUpperNodes != nLowerNodes || nUpperElems != nLowerElems) {
		std::cout << "***Error: nUpperNodes != nLowerNodes or nUpperElems != nLowerElems\n";
		return nullptr;
	}

	double lower_nodes_max_z = -std::numeric_limits<double>::max();
	double lower_nodes_min_z = std::numeric_limits<double>::max();
	for (auto itr=lower_nodes.begin(); itr!=lower_nodes.end(); ++itr) {
		lower_nodes_max_z = std::max(lower_nodes_max_z, (*itr)->getCoords()[2]);
		lower_nodes_min_z = std::min(lower_nodes_min_z, (*itr)->getCoords()[2]);
	}
    double upper_nodes_max_z = -std::numeric_limits<double>::max();
    double upper_nodes_min_z = std::numeric_limits<double>::max();
    for (auto itr=lower_nodes.begin(); itr!=lower_nodes.end(); ++itr) {
        lower_nodes_max_z = std::max(lower_nodes_max_z, (*itr)->getCoords()[2]);
        lower_nodes_min_z = std::min(lower_nodes_min_z, (*itr)->getCoords()[2]);
    }
    for (auto itr=upper_nodes.begin(); itr!=upper_nodes.end(); ++itr) {
        upper_nodes_max_z = std::max(upper_nodes_max_z, (*itr)->getCoords()[2]);
        upper_nodes_min_z = std::min(upper_nodes_min_z, (*itr)->getCoords()[2]);
    }
	double min_z = lower_nodes_min_z - extZ;
	double max_z = lower_nodes_max_z + dz + extZ;

//#lower
//	size_t n_new_nodes = nUpperNodes + (nHalfLayers * nUpperNodes);
//	size_t n_new_eles = nUpperElems * nHalfLayers; //( 1 + nHalfLayers);
//#lower + connect
//	size_t n_new_nodes = nUpperNodes*2 + (nHalfLayers * nUpperNodes);
//	size_t n_new_eles = nUpperElems * (1 + nHalfLayers);
//#all
	size_t n_new_nodes = nUpperNodes * (2 + 2 * nHalfLayers);
	size_t n_new_eles = nUpperElems * ( 1 + 2 * nHalfLayers);

    std::cout << "\n======================================\n";
    std::cout << "#statistics:\n";
    std::cout << "lower surface z min, max: " << lower_nodes_min_z << ", " << lower_nodes_max_z << "\n";
    std::cout << "upper surface z min, max: " << upper_nodes_min_z << ", " << upper_nodes_max_z << "\n";
	std::cout << "\n======================================\n";
    std::cout << "#given parameters:\n";
	std::cout << "Z extension in half: " << extZ << "\n";
	std::cout << "delta Z: " << dz << "\n";

	std::vector<MeshLib::Node*> new_nodes(n_new_nodes);
	std::vector<MeshLib::Element*> new_elems(n_new_eles);
	size_t node_offset = 0;
	size_t elem_offset = 0;
	double z_offset(0.0);
	unsigned mat_id = 0;
	bool mat_id_each_layer = false;


	addHalfLayers(lower_nodes, lower_elems, nHalfLayers, z_offset, dz, min_z, mat_id, mat_id_each_layer, new_nodes, new_elems, node_offset, elem_offset);
	unsigned mat_id_connect = mat_id + 1;
	connect(lower_nodes, lower_elems, upper_nodes, upper_elems, dz, mat_id_connect, new_nodes, new_elems, node_offset, elem_offset);
	z_offset = -dz;
	//if (mat_id_each_layer)
	    mat_id = mat_id + 2;
	addHalfLayers(upper_nodes, upper_elems, nHalfLayers, z_offset, -dz, max_z, mat_id, mat_id_each_layer, new_nodes, new_elems, node_offset, elem_offset);

    std::cout << "\n======================================\n";
    std::cout << "#result:\n";
    std::cout << "nr. of half layers: " << nHalfLayers << "\n";
    std::cout << "volume z min, max: " << min_z << ", " << max_z << "\n";
    std::cout << "nr. of surface nodes and elements: " << nUpperNodes << ", " << nUpperElems << "\n";
    std::cout << "nr. of expected new nodes and elements: " << n_new_nodes << ", " << n_new_eles << "\n";

	return new MeshLib::Mesh("SubsurfaceMesh", new_nodes, new_elems);
}

MeshLib::Mesh* CreateMinLayers(const MeshLib::Mesh* uppesr_mesh, const MeshLib::Mesh* lower_mesh, std::size_t nHalfLayers, double dz)
{
    const size_t nUpperNodes = uppesr_mesh->getNNodes();
    const size_t nUpperElems = uppesr_mesh->getNElements();
    const std::vector<MeshLib::Node*> upper_nodes = uppesr_mesh->getNodes();
    const std::vector<MeshLib::Element*> upper_elems = uppesr_mesh->getElements();
    const size_t nLowerNodes = lower_mesh->getNNodes();
    const size_t nLowerElems = lower_mesh->getNElements();
    const std::vector<MeshLib::Node*> lower_nodes = lower_mesh->getNodes();
    const std::vector<MeshLib::Element*> lower_elems = lower_mesh->getElements();

    if (nUpperNodes != nLowerNodes || nUpperElems != nLowerElems) {
        ERR("***Error: nUpperNodes != nLowerNodes or nUpperElems != nLowerElems");
        return nullptr;
    }

    double lower_nodes_max_z = -std::numeric_limits<double>::max();
    double lower_nodes_min_z = std::numeric_limits<double>::max();
    for (auto itr=lower_nodes.begin(); itr!=lower_nodes.end(); ++itr) {
        lower_nodes_max_z = std::max(lower_nodes_max_z, (*itr)->getCoords()[2]);
        lower_nodes_min_z = std::min(lower_nodes_min_z, (*itr)->getCoords()[2]);
    }
    double upper_nodes_max_z = -std::numeric_limits<double>::max();
    double upper_nodes_min_z = std::numeric_limits<double>::max();
    for (auto itr=lower_nodes.begin(); itr!=lower_nodes.end(); ++itr) {
        lower_nodes_max_z = std::max(lower_nodes_max_z, (*itr)->getCoords()[2]);
        lower_nodes_min_z = std::min(lower_nodes_min_z, (*itr)->getCoords()[2]);
    }
    for (auto itr=upper_nodes.begin(); itr!=upper_nodes.end(); ++itr) {
        upper_nodes_max_z = std::max(upper_nodes_max_z, (*itr)->getCoords()[2]);
        upper_nodes_min_z = std::min(upper_nodes_min_z, (*itr)->getCoords()[2]);
    }
    double diff_z_min = std::numeric_limits<double>::max();
    double diff_z_max = -std::numeric_limits<double>::max();
    for (std::size_t i=0; i<nUpperNodes; i++) {
        diff_z_min = std::min(diff_z_min, upper_nodes[i]->getCoords()[2]-lower_nodes[i]->getCoords()[2]);
        diff_z_max = std::max(diff_z_max, upper_nodes[i]->getCoords()[2]-lower_nodes[i]->getCoords()[2]);
    }

//#lower
//  size_t n_new_nodes = nUpperNodes + (nHalfLayers * nUpperNodes);
//  size_t n_new_eles = nUpperElems * nHalfLayers; //( 1 + nHalfLayers);
//#lower + connect
//  size_t n_new_nodes = nUpperNodes*2 + (nHalfLayers * nUpperNodes);
//  size_t n_new_eles = nUpperElems * (1 + nHalfLayers);
//#all
#ifdef BOTH_SIDES
    size_t n_new_nodes = nUpperNodes * (2 + nHalfLayers + 1);
    size_t n_new_eles = nUpperElems * ( 1 + nHalfLayers + 1);
#else
    size_t n_new_nodes = nUpperNodes * (1 + nHalfLayers + 1);
    size_t n_new_eles = nUpperElems * ( nHalfLayers + 1);
#endif

    std::cout << "\n======================================\n";
    std::cout << "#statistics:\n";
    std::cout << "lower surface z min, max: " << lower_nodes_min_z << ", " << lower_nodes_max_z << "\n";
    std::cout << "upper surface z min, max: " << upper_nodes_min_z << ", " << upper_nodes_max_z << "\n";
    std::cout << "difference min, max: " << diff_z_min << ", " << diff_z_max << "\n";
    std::cout << "\n======================================\n";
    std::cout << "#given parameters:\n";
    std::cout << "nr. of half layers: " << nHalfLayers << "\n";
    std::cout << "delta Z: " << dz << "\n";

    double void_dz =  std::abs(diff_z_min)*1.1;
//    if (diff_z_min < .0 && dz < std::abs(diff_z_min)*1.1) {
//        dz = std::abs(diff_z_min)*1.1;
//        std::cout << "-> delta Z is adjusted to : " << dz << "\n";
//    }

    std::vector<MeshLib::Node*> new_nodes(n_new_nodes);
    std::vector<MeshLib::Element*> new_elems(n_new_eles);
    size_t node_offset = 0;
    size_t elem_offset = 0;
    double z_offset(0.0);
    const bool mat_id_each_layer = false;


    // only one layer for bottom
    unsigned mat_id = 0;
#ifdef BOTH_SIDES
    INFO("create bottom layers...");
    extrudeLayer(lower_nodes, lower_elems, 1, z_offset, dz, mat_id, mat_id_each_layer, new_nodes, new_elems, node_offset, elem_offset);
    mat_id++;
 #endif
    // add interface
    INFO("create interface layers...");
#ifdef BOTH_SIDES
    unsigned mat_id_connect = 2;
    connect(lower_nodes, lower_elems, upper_nodes, upper_elems, dz, mat_id_connect, new_nodes, new_elems, node_offset, elem_offset);
#else
    unsigned mat_id_connect = 1;
    createVoidLayer(lower_nodes, lower_elems, upper_nodes, upper_elems, void_dz, mat_id_connect, new_nodes, new_elems, node_offset, elem_offset);
#endif
    z_offset = -void_dz;
    // add interface
    INFO("create upper layers...");
    extrudeLayer(upper_nodes, upper_elems, nHalfLayers, z_offset, -dz, mat_id, mat_id_each_layer, new_nodes, new_elems, node_offset, elem_offset);
    mat_id++;

    std::cout << "\n======================================\n";
    std::cout << "#result:\n";
    std::cout << "nr. of half layers: " << nHalfLayers << "\n";
    std::cout << "nr. of surface nodes and elements: " << nUpperNodes << ", " << nUpperElems << "\n";
    std::cout << "nr. of expected new nodes and elements: " << n_new_nodes << ", " << n_new_eles << "\n";

    return new MeshLib::Mesh("SubsurfaceMesh", new_nodes, new_elems);
}


int main(int argc, char* argv[])
{
    LOGOG_INITIALIZE();
    logog::Cout* logogCout = new logog::Cout;
    BaseLib::LogogSimpleFormatter* formatter = new BaseLib::LogogSimpleFormatter;
    logogCout->SetFormatter(*formatter);

    try{
    TCLAP::CmdLine cmd("Create a 3D layer mesh from two surface meshes (for DECOVALEX Task C1).", ' ', "0.1");
    TCLAP::ValueArg<std::string> arg_mesh_upper("u", "upper-mesh-file",
                                         "upper mesh file", true,
                                         "", "file name of upper mesh");
    cmd.add(arg_mesh_upper);
    TCLAP::ValueArg<std::string> arg_mesh_lower("l", "lower-mesh-file",
                                         "lower mesh file", true,
                                         "", "file name of lower mesh");
    cmd.add(arg_mesh_lower);
    TCLAP::ValueArg<std::string> arg_mesh_out("o", "mesh-output-file",
                                          "the name of the file the mesh will be written to", true,
                                          "", "file name of output mesh");
    cmd.add(arg_mesh_out);
    TCLAP::ValueArg<double> arg_halfZExt("e", "half-z-extension", "half Z extension", false,
                                          1, "half Z extension");
    cmd.add(arg_halfZExt);
    TCLAP::ValueArg<unsigned> arg_nLayers("n", "nr-layers", "the number of layers", false,
                                          1, "the number of layers");
    cmd.add(arg_nLayers);
    TCLAP::ValueArg<double> arg_dz("d", "dz", "delta z", true, 0.1, "delta z");
    cmd.add(arg_dz);

    TCLAP::SwitchArg arg_minimum("m", "minimum", "only connected", false);
    cmd.add(arg_minimum);

    TCLAP::ValueArg<std::string> arg_direct_file("c", "condition-direct-file",
                                          "the name of the file direct conditions will be written to", false,
                                          "", "file name of output direct");
    cmd.add(arg_direct_file);

    cmd.parse(argc, argv);


	std::string upper_msh_name(arg_mesh_upper.getValue()), lower_msh_name(arg_mesh_lower.getValue()), new_mesh_name(arg_mesh_out.getValue());

	MeshLib::Mesh* upper_mesh (FileIO::readMeshFromFile(upper_msh_name));
	MeshLib::Mesh* lower_mesh (FileIO::readMeshFromFile(lower_msh_name));

	double extZ = arg_halfZExt.getValue();
	double dz = arg_dz.getValue();
	unsigned nHalfLayers = arg_nLayers.getValue();


    MeshLib::Mesh* mesh_new = nullptr;
    if (arg_minimum.isSet())
        mesh_new = CreateMinLayers(upper_mesh, lower_mesh, nHalfLayers, dz);
    else
        mesh_new = CreateLayers(upper_mesh, lower_mesh, extZ, dz);

	if (mesh_new) {
		std::cout << "\n======================================\n";
		std::cout << "Generated mesh\n";
		std::cout << "nr. of nodes: " << mesh_new->getNNodes() << "\n";
		std::cout << "nr. of elements: " << mesh_new->getNElements() << "\n";
		std::cout << "\n======================================\n";

		FileIO::VtuInterface meshIO(mesh_new);
//		meshIO.setPrecision(9);
		meshIO.writeToFile(new_mesh_name);

		if (arg_direct_file.isSet()) {
			std::cout << "output direct BCs into a file:" << arg_direct_file.getValue() << std::endl;
			std::ofstream of(arg_direct_file.getValue().c_str());
			if (of) {
				for (size_t i=0; i<lower_mesh->getNNodes(); i++) {
					of << i << " 0\n";
				}
			}
			of.close();
		}
	}

	} catch (TCLAP::ArgException &e) {
	    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

    delete formatter;
    delete logogCout;
    LOGOG_SHUTDOWN();

    return 0;
}


