/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "ComparisonTools.h"

#include <vector>
#include <cmath>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

using namespace boost::property_tree;

namespace
{

struct MeshInfo
{
    size_t n_nodes;
    size_t n_elements;

    bool operator ==(const MeshInfo &msh) {
        if (n_nodes != msh.n_nodes) return false;
        if (n_elements != msh.n_elements) return false;
        return true;
    }
    bool operator !=(const MeshInfo &msh) {
        return !((*this)==msh);
    }
};

void getMeshInfo(const ptree &doc, MeshInfo &msh)
{
    auto &e_piece = doc.get_child("VTKFile.UnstructuredGrid.Piece");
    auto &attributes = e_piece.get_child("<xmlattr>");
    msh.n_nodes = attributes.get<std::size_t>("NumberOfPoints");
    msh.n_elements = attributes.get<std::size_t>("NumberOfCells");
}

void getListofDataArray(const ptree& e_ptdata, const size_t n_pt, std::vector<std::vector<double> > &vec_pointdata, std::vector<std::string> &vec_pointdata_name)
{
    size_t i_array = 0;
    auto range = e_ptdata.equal_range("DataArray");
    for (auto it=range.first; it!=range.second; ++it)
    {
        auto &attributes = it->second.get_child("<xmlattr>");
        vec_pointdata_name.push_back(attributes.get<std::string>("Name", ""));
        unsigned n_comp = attributes.get<unsigned>("NumberOfComponents", 1u);
        size_t n_total_len = n_comp * n_pt;
        vec_pointdata.push_back(std::vector<double>());
        vec_pointdata[i_array].resize(n_total_len);
        std::vector<double> &data = vec_pointdata[i_array];

        std::stringstream ss(std::stringstream::in | std::stringstream::out);
        ss << it->second.get_value<std::string>();

        for (size_t i=0; i<n_total_len; i++) {
            ss >> data[i];
        }

        i_array++;
    }
}

void getPointData(ptree &doc, size_t n_pt, std::vector<std::vector<double> > &vec_pointdata, std::vector<std::string> &vec_pointdata_name)
{
    auto& e_ptdata = doc.get_child("VTKFile.UnstructuredGrid.Piece.PointData");
    getListofDataArray(e_ptdata, n_pt, vec_pointdata, vec_pointdata_name);
}


void getCellData(ptree &doc, size_t n_pt, std::vector<std::vector<double> > &vec_pointdata, std::vector<std::string> &vec_pointdata_name)
{
    auto& e_celldata = doc.get_child("VTKFile.UnstructuredGrid.Piece.CellData");
    getListofDataArray(e_celldata, n_pt, vec_pointdata, vec_pointdata_name);
}

bool compareRealVector(const std::vector<double> &ref, const std::vector<double> &test, const double epsilon)
{
    const size_t n = ref.size();
    size_t cnt_fail = 0;
    double max_val = .0;
    for (size_t i=0; i<n; i++) {
        max_val = std::max(max_val, std::abs(ref[i]));
    }
    double rel_error = epsilon;
    if (max_val != .0)
        rel_error *= max_val;
    std::cout << "\t rel. tol. = " << epsilon << ", abs. tol. =" << rel_error << std::endl;
    for (size_t i=0; i<n; i++) {
        double diff = ref[i] - test[i];
        if (std::abs(diff) > rel_error) {
            cnt_fail++;
            std::cout << "*** DIFF in " << i << ": Ref=" << ref[i] << ", New=" << test[i] << std::endl;
        }
    }

    return cnt_fail==0;
}

} //end namespace

bool compareAciiVTUFiles(const std::string &str_ref_file, const std::string &str_test_file, const double epsilon)
{
    // check if files exist
    if (!BaseLib::IsFileExisting(str_ref_file)) {
        std::cout << "Reference file not found: " << str_ref_file << std::endl;
        return false;
    }
    if (!BaseLib::IsFileExisting(str_test_file)) {
        std::cout << "New result file not found: " << str_test_file << std::endl;
        return false;
    }

    // open both XML files
    ptree doc_ref;
    read_xml(str_ref_file, doc_ref);
    ptree doc_test;
    read_xml(str_test_file, doc_test);

    // compare contents

    // check mesh
    MeshInfo msh_ref, msh_test;
    getMeshInfo(doc_ref, msh_ref);
    getMeshInfo(doc_test, msh_test);
    if (msh_ref != msh_test) {
        std::cout << "*** DIFF:" << std::endl;
        std::cout << "Ref: " << msh_ref.n_nodes << " nodes, " << msh_ref.n_elements << " elements" << std::endl;
        std::cout << "New: " << msh_test.n_nodes << " nodes, " << msh_test.n_elements << " elements" << std::endl;
        return false;
    }

    // check point data
    {
        size_t n_pt = msh_ref.n_nodes;
        std::vector<std::vector<double> > vec_pointdata_ref;
        std::vector<std::vector<double> > vec_pointdata_test;
        std::vector<std::string> vec_pointdata_name_ref, vec_pointdata_name_test;
        getPointData(doc_ref, n_pt, vec_pointdata_ref, vec_pointdata_name_ref);
        getPointData(doc_test, n_pt, vec_pointdata_test, vec_pointdata_name_test);
        if (vec_pointdata_ref.size() != vec_pointdata_test.size()) {
            std::cout << "*** DIFF:" << std::endl;
            std::cout << "Ref: " << vec_pointdata_ref.size() << " point data" << std::endl;
            std::cout << "New: " << vec_pointdata_test.size() << " point data" << std::endl;
            return false;
        }
        for (size_t i=0; i<vec_pointdata_ref.size(); i++) {
            std::cout << "-->PointData: " << vec_pointdata_name_ref[i] << std::endl;
            if (!compareRealVector(vec_pointdata_ref[i], vec_pointdata_test[i], epsilon))
                return false;
        }
    }

    // check cell data
    size_t n_cell = msh_ref.n_elements;
    std::vector<std::vector<double> > vec_celldata_ref;
    std::vector<std::vector<double> > vec_celldata_test;
    std::vector<std::string> vec_celldata_name_ref, vec_celldata_name_test;
    getCellData(doc_ref, n_cell, vec_celldata_ref, vec_celldata_name_ref);
    getCellData(doc_test, n_cell, vec_celldata_test, vec_celldata_name_test);
    if (vec_celldata_ref.size() != vec_celldata_test.size()) {
        std::cout << "*** DIFF:" << std::endl;
        std::cout << "Ref: " << vec_celldata_ref.size() << " cell data" << std::endl;
        std::cout << "New: " << vec_celldata_test.size() << " cell data" << std::endl;
        return false;
    }
    for (size_t i=0; i<vec_celldata_ref.size(); i++) {
        if (vec_celldata_name_ref[i] == "MatGroup") continue;
        std::cout << "-->CellData: " << vec_celldata_name_ref[i] << std::endl;
        if (!compareRealVector(vec_celldata_ref[i], vec_celldata_test[i], epsilon)) {
            return false;
        }
    }
    return true;
}
