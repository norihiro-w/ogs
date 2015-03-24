/**
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <array>
#include <string>
#include <valarray>

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"

#include "BaseLib/BuildInfo.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"

#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"

#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"

#include "FileIO/readMeshFromFile.h"
#include "FileIO/writeMeshToFile.h"

MeshLib::Mesh* addInletAndOutlet(MeshLib::Mesh* org_mesh)
{
	// define lines to define inlet and outlet boundaries
	const double L = 0.0895;
	const double W = 0.05;
	std::vector<GeoLib::Point*> vec_pt1;
	vec_pt1.push_back(new GeoLib::Point(0, 0, 0));
	vec_pt1.push_back(new GeoLib::Point(0, W, 0));
	std::vector<GeoLib::Point*> vec_pt2;
	vec_pt2.push_back(new GeoLib::Point(L, 0, 0));
	vec_pt2.push_back(new GeoLib::Point(L, W, 0));
	GeoLib::Polyline *ply1 = new GeoLib::Polyline(vec_pt1);
	ply1->addPoint(0);
	ply1->addPoint(1);
	GeoLib::Polyline *ply2 = new GeoLib::Polyline(vec_pt2);
	ply2->addPoint(0);
	ply2->addPoint(1);

	// pick up nodes at inlet and outlet boundaries
	MeshGeoToolsLib::MeshNodesAlongPolyline nodeSearch1(*org_mesh, *ply1, 1e-6);
	std::vector<std::size_t> nodeIDs_left = nodeSearch1.getNodeIDs();
	MeshGeoToolsLib::MeshNodesAlongPolyline nodeSearch2(*org_mesh, *ply2, 1e-6);
	std::vector<std::size_t> nodeIDs_right = nodeSearch2.getNodeIDs();

	std::cout << "-> " << nodeIDs_left.size() << " nodes found at the left boundary" << std::endl;
	std::cout << "-> " << nodeIDs_right.size() << " nodes found at the right boundary" << std::endl;

	// copy mesh
	std::vector<MeshLib::Node*> new_nodes(MeshLib::copyNodeVector(org_mesh->getNodes()));
	std::vector<MeshLib::Element*> new_eles(MeshLib::copyElementVector(org_mesh->getElements(), new_nodes));
	std::size_t node_counter = new_nodes.size();
	std::size_t ele_counter = new_eles.size();

	// generate nodes for extra elements (4 layers)
	//const auto n_org_nodes = new_nodes.size();
	const unsigned n_layers = 5u;
	const double d = 1e-3;
	std::cout << "-> create new elements with " << n_layers << " layers and increment length of " << d  << std::endl;
	// inlet
	const std::valarray<double> n_vec1({-1, 0, 0});
	std::vector<MeshLib::Node*> new_nodes_inlet;
	for (auto node_id : nodeIDs_left) {
		MeshLib::Node& org_node = *new_nodes[node_id];
		for (unsigned il=n_layers; il>0; il--) {
			auto dv = n_vec1 * d*il;
			MeshLib::Node* new_node = new MeshLib::Node(org_node[0]+dv[0],org_node[1]+dv[1],org_node[2]+dv[2], node_counter++);
			new_nodes.push_back(new_node);
			new_nodes_inlet.push_back(new_node);
		}
		new_nodes_inlet.push_back(&org_node);
	}

	// outlet
	const std::valarray<double> n_vec2({1, 0, 0});
	std::vector<MeshLib::Node*> new_nodes_outlet;
	for (auto node_id : nodeIDs_right) {
		MeshLib::Node& org_node = *new_nodes[node_id];
		new_nodes_outlet.push_back(&org_node);
		for (unsigned il=0; il<n_layers; il++) {
			auto dv = n_vec2 * d*(il+1);
			MeshLib::Node* new_node = new MeshLib::Node(org_node[0]+dv[0],org_node[1]+dv[1],org_node[2]+dv[2], node_counter++);
			new_nodes.push_back(new_node);
			new_nodes_outlet.push_back(new_node);
		}
	}

	// construct new elements
	// inlet
	{
		const std::size_t n_y_cells = nodeIDs_left.size() - 1;
		const std::size_t n_x_cells = n_layers;
		const std::size_t n_x_nodes = n_x_cells + 1;
		for (std::size_t j = 0; j < n_y_cells; j++)
		{
			const std::size_t offset_y1 = j * n_x_nodes;
			const std::size_t offset_y2 = (j + 1) * n_x_nodes;
			for (std::size_t k = 0; k < n_x_cells; k++)
			{
				std::array<MeshLib::Node*, 3> element1_nodes;
				element1_nodes[0] = new_nodes_inlet[offset_y1 + k];
				element1_nodes[1] = new_nodes_inlet[offset_y2 + k + 1];
				element1_nodes[2] = new_nodes_inlet[offset_y2 + k];
				new_eles.push_back (new MeshLib::Tri(element1_nodes, 1, ele_counter++));
				std::array<MeshLib::Node*, 3> element2_nodes;
				element2_nodes[0] = new_nodes_inlet[offset_y1 + k];
				element2_nodes[1] = new_nodes_inlet[offset_y1 + k + 1];
				element2_nodes[2] = new_nodes_inlet[offset_y2 + k + 1];
				new_eles.push_back (new MeshLib::Tri(element2_nodes, 1, ele_counter++));
			}
		}
	}

	// outlet
	{
		const std::size_t n_y_cells = nodeIDs_right.size() - 1;
		const std::size_t n_x_cells = n_layers;
		const std::size_t n_x_nodes = n_x_cells + 1;
		for (std::size_t j = 0; j < n_y_cells; j++)
		{
			const std::size_t offset_y1 = j * n_x_nodes;
			const std::size_t offset_y2 = (j + 1) * n_x_nodes;
			for (std::size_t k = 0; k < n_x_cells; k++)
			{
				std::array<MeshLib::Node*, 3> element1_nodes;
				element1_nodes[0] = new_nodes_outlet[offset_y1 + k];
				element1_nodes[1] = new_nodes_outlet[offset_y2 + k + 1];
				element1_nodes[2] = new_nodes_outlet[offset_y2 + k];
				new_eles.push_back (new MeshLib::Tri(element1_nodes, 1, ele_counter++));
				std::array<MeshLib::Node*, 3> element2_nodes;
				element2_nodes[0] = new_nodes_outlet[offset_y1 + k];
				element2_nodes[1] = new_nodes_outlet[offset_y1 + k + 1];
				element2_nodes[2] = new_nodes_outlet[offset_y2 + k + 1];
				new_eles.push_back (new MeshLib::Tri(element2_nodes, 1, ele_counter++));
			}
		}
	}


	// construct new mesh
	return new MeshLib::Mesh(org_mesh->getName() + "_new", new_nodes, new_eles);
}

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Task C1 tool - add extra elements for inlet and outlet", ' ', BaseLib::BuildInfo::git_describe);
	TCLAP::ValueArg<std::string> input_arg("i", "input-mesh-file","input mesh file",true,"","string");
	cmd.add( input_arg );
	TCLAP::ValueArg<std::string> output_arg("o", "output-mesh-file","output mesh file",true,"","string");
	cmd.add( output_arg );
	cmd.parse( argc, argv );

	// read a mesh file
	MeshLib::Mesh* org_mesh (FileIO::readMeshFromFile(input_arg.getValue()));
	if (!org_mesh)
		return 0;
	INFO("Mesh read: %d nodes, %d elements.", org_mesh->getNNodes(), org_mesh->getNElements());

	// revise the mesh
	MeshLib::Mesh* new_mesh = addInletAndOutlet(org_mesh);

	// write into a file
	if (new_mesh) {
		INFO("New mesh: %d nodes, %d elements.", new_mesh->getNNodes(), new_mesh->getNElements());
		FileIO::writeMeshToFile(*new_mesh, output_arg.getValue());
	}

	delete org_mesh;
	delete new_mesh;

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
