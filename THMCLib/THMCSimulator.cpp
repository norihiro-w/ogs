/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "THMCSimulator.h"

#include <algorithm>
#include <iostream>

// external library
#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <logog/include/logog.hpp>
#include <tclap/CmdLine.h>
#ifdef USE_LIS
#include <lis.h>
#endif
#ifdef USE_PETSC
#include <petscksp.h>
#endif

// internal library
#include "BaseLib/CodingTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/LogogCustomCout.h"
#include "BaseLib/TemplateLogogFormatterSuppressedGCC.h"
#include "BaseLib/MPITools.h"

// NumLib
#include "NumLib/TransientCoupling/TransientCouplingStructureBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"

#include "FileIO/Legacy/Ogs5FemIO.h"

// this module
#include "SimulationInfo.h"
#include "Ogs6FemData.h"
#include "Ogs5ToOgs6.h"
#include "TimeSteppingControllerWithOutput.h"

namespace ogs6
{

//##############################################################################
// Variables
//##############################################################################
static logog::Target *logogCout;
static logog::LogFile *logog_file;
static logog::Formatter *custom_format;
static bool isOgsInitCalled = false;
static bool isOgsExitCalled = false;

//##############################################################################
// Functions
//##############################################################################
void ogsInit(int argc, char* argv[])
{
    if (isOgsInitCalled) return;
    isOgsInitCalled = true;

#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif

    LOGOG_INITIALIZE();
    custom_format = new BaseLib::TemplateLogogFormatterSuppressedGCC<TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG>(); //BaseLib::LogogSimpleFormatter();
    logogCout = new BaseLib::LogogCustomCout(); //new logog::Cout();
    logogCout->SetFormatter(*custom_format);
    logog_file = NULL;

#ifdef USE_LIS
    lis_initialize(&argc, &argv);
#endif
#ifdef USE_PETSC
    char petsc_help[] = "Using PETSc package\n";
    PetscInitialize(&argc, &argv,(char *)0,petsc_help);
#endif
}

void ogsExit()
{
    if (isOgsExitCalled) return;
    isOgsExitCalled = true;

#ifdef USE_LIS
    lis_finalize();
#endif
#ifdef USE_PETSC
    PetscFinalize();
#endif

    INFO("exit ogs6.");
    delete custom_format;
    delete logogCout;
    delete logog_file;
    LOGOG_SHUTDOWN();

#ifdef USE_MPI
    MPI_Finalize();
#endif
}

//##############################################################################
// THMCSimulator
//##############################################################################
THMCSimulator::THMCSimulator(int argc, char* argv[])
: _sim_info(NULL), _cpl_system(NULL)
{
    try {
        // Command line parser
        TCLAP::CmdLine cmd("ogs6", ' ', "0.1");
        cmd.ignoreUnmatched(true);
        TCLAP::ValueArg<std::string> input_arg("i", "input", "input file", false, "", "string");
        cmd.add( input_arg );
        TCLAP::ValueArg<std::string> output_dir_arg("o", "output", "output directory", false, "", "string");
        cmd.add( output_dir_arg );
        TCLAP::ValueArg<std::string> logfile_arg("l", "log", "log file", false, "", "string");
        cmd.add( logfile_arg );
        TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "number");
        cmd.add( verbosity_arg );
        TCLAP::SwitchArg pcs_arg("m", "modules", "list available modules", false);
        cmd.add( pcs_arg );
        TCLAP::ValueArg<std::string> alglib_arg("", "ls", "linear algebraic library", false, "eigen", "eigen|eigen_lis|lis|petsc");
        cmd.add( alglib_arg );
        cmd.parse( argc, argv ); // process can exit in this function

        // initialize
        ogsInit(argc, argv);

        // get parsed data
        // log file
        if (! logfile_arg.getValue().empty()) {
            if (!logog_file) delete logog_file;
            std::string log_file = logfile_arg.getValue();
            BaseLib::truncateFile(log_file); // do this not to append log into an existing file
            logog_file = new logog::LogFile(log_file.c_str());
            logog_file->SetFormatter( *custom_format );
        }

        SimulationInfo::outputHeader();
        // list modules
        const unsigned flag_list_modules (pcs_arg.getValue());
        if (flag_list_modules!=0) {
            ProcessBuilder::getInstance()->output();
        }

        const bool is_input_file_given = !input_arg.getValue().empty();
        if (is_input_file_given) {
            INFO("->Parsing input arguments");
            INFO("* project path     : %s", input_arg.getValue().c_str());
            if (! logfile_arg.getValue().empty()) {
                INFO("* log file path    : %s", logfile_arg.getValue().c_str());
            }

            // data output directory
            std::string output_dir_path = "";
            if (! output_dir_arg.getValue().empty()) {
                output_dir_path = output_dir_arg.getValue();
            }
            INFO("* output directory : %s", output_dir_path.c_str());

            if (! input_arg.getValue().empty()) {
                const std::string proj_path = input_arg.getValue();
                if (checkInputFiles(proj_path)) {
                    _sim_info = new SimulationInfo(proj_path, output_dir_path);
                } else {
                    ERR("***Error: Cannot find a project - %s", proj_path.c_str());
                }
            }
        }

        if (alglib_arg.isSet()) {
            INFO("* linear solver library : %s", alglib_arg.getValue().c_str());
            if (alglib_arg.getValue()=="eigen")
                THMCLib::Ogs6FemData::getInstance()->linalg_type = MathLib::LinAlgLibType::Eigen;
#ifdef USE_LIS
            else if (alglib_arg.getValue()=="eigen_lis")
                THMCLib::Ogs6FemData::getInstance()->linalg_type = MathLib::LinAlgLibType::EigenLis;
            else if (alglib_arg.getValue()=="lis")
                THMCLib::Ogs6FemData::getInstance()->linalg_type = MathLib::LinAlgLibType::Lis;
#endif
#ifdef USE_PETSC
            else if (alglib_arg.getValue()=="petsc")
                THMCLib::Ogs6FemData::getInstance()->linalg_type = MathLib::LinAlgLibType::PETSc;
#endif
            else
                ERR("Unsupported linear algebraic library %s", alglib_arg.getValue().c_str());
        }


    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

}

THMCSimulator::~THMCSimulator()
{
    THMCLib::Ogs6FemData::getInstance()->initialize();
    delete _sim_info;
    delete _cpl_system;
    ogsExit();
}

bool THMCSimulator::checkInputFiles(const std::string& proj_path)
{
    // meanwhile OGS5 files are default
    if(!BaseLib::IsFileExisting(proj_path+".pcs"))
    {
        ERR("Cannot find a PCS file - %s.pcs", proj_path.c_str());
        return false;
    }
//    if(!BaseLib::IsFileExisting(proj_path+ ".pro"))
//    {
//        ERR("Cannot find a property file - %s.pro", proj_path.c_str());
//        return false;
//    }

    return true;
}

int THMCSimulator::execute()
{
    if (!_sim_info) return 0;

    auto* ogs6fem = THMCLib::Ogs6FemData::getInstance();
    const std::string proj_path = _sim_info->getProjectPath();
    ogs6fem->project_name = _sim_info->getProjectName();
    ogs6fem->output_dir = _sim_info->getOutputDirPath();

    //-------------------------------------------------------------------------
    // Read files
    //-------------------------------------------------------------------------
    INFO("->Reading input files...");
    // coupling
    using namespace boost::property_tree;
    ptree opRoot;
    auto &opOgs = opRoot.put_child("ogs", ptree());
    // ogs5fem
    if (!BaseLib::IsFileExisting(proj_path+".pcs")) {
    	return 0;
    }
    INFO("->Reading OGS5 input files...");
    ogs5::Ogs5FemData ogs5femdata;
    INFO("-------------------------------------------------");
    ogs5::Ogs5FemIO::read(proj_path, ogs5femdata);
    INFO("-------------------------------------------------");
    if (!Ogs5ToOgs6::convert(ogs5femdata, *ogs6fem, opOgs)) {
        ERR("***Error: Failure during conversion of ogs5 to ogs6.");
        return 0;
    }

    // output converted setting
    BaseLib::MPIEnvironment mpi;
    if (mpi.root()) {
        std::string str_conversion_logfile = (ogs6fem->output_dir.empty() ? "." : ogs6fem->output_dir) + "/converted_setting.log";
#if BOOST_VERSION >= 105600
        auto settings = boost::property_tree::xml_writer_make_settings<std::string> ('\t', 1);
#else
        xml_writer_settings<char> settings('\t', 1);
#endif
        write_xml(str_conversion_logfile, opRoot, std::locale(), settings);
    }


	//-------------------------------------------------------------------------
    // Setup simulation
    //-------------------------------------------------------------------------
    // construct mesh
    INFO("->Constructing meshes... %d mesh loaded", ogs6fem->list_mesh.size());
    for (size_t i=0; i<ogs6fem->list_mesh.size(); i++) {
        MeshLib::Mesh* msh = ogs6fem->list_mesh[i];
        //msh->constructGeometricProperty();
        INFOa("->mesh id %d: dim=%d, nodes=%d, elements=%d", i, msh->getDimension(), msh->getNNodes(), msh->getNElements());
    }

    // construct element data for FE calculations
    ogs6fem->fe_mesh_data.resize(ogs6fem->list_mesh.size());
//    ogs6fem->fe_ele_data.resize(ogs6fem->list_mesh.size());
//    for (size_t i=0; i<ogs6fem->list_mesh.size(); i++) {
//        MeshLib::Mesh* msh = ogs6fem->list_mesh[i];
//        ogs6fem->fe_ele_data[i].resize(msh->getNElements());
//    }

    // construct coupling system
    INFO("->Generating coupling system...");
    typedef class NumLib::TemplateCouplingStrucutreBuilder
        <
        NumLib::ITransientCoupledSystem,
		THMCLib::Process,
        NumLib::AsyncPartitionedSystem,
        NumLib::TransientPartitionedAlgorithmFactory
        > CoupledProcessStrucutreBuilder;

    CoupledProcessStrucutreBuilder cpl_builder;
    delete _cpl_system;
    _cpl_system = cpl_builder.build(&opOgs, *GeoProcessBuilder::getInstance());
    std::vector<std::string> &list_mono_system_name = cpl_builder.getListOfMonolithicSystemName();
    if (list_mono_system_name.size()==0) {
        ERR("***Error: no active process is selected.");
        return 0;
    }
    if (!_cpl_system->check()) {
        ERR("***Error while checking coupled system");
        return 0;
    }

    // list up monolithic processes
    INFO("->Initializing all processes...");
    std::vector<THMCLib::Process*> &list_mono_system = cpl_builder.getListOfMonolithicSystem();
    for (size_t i=0; i<list_mono_system.size(); i++) {
        std::string &pcs_name = list_mono_system_name[i];
        THMCLib::Process* pcs = list_mono_system[i];
        if (pcs->getProcessName().empty())
            pcs->setProcessName(pcs->getProcessType());
        INFO("PCS %d: name=%s, type=%s (IN=%d, OUT=%d)", i, pcs_name.c_str(), pcs->getProcessType().c_str(), pcs->getNumberOfInputParameters(), pcs->getNumberOfOutputParameters());
        for (size_t j=0; j<pcs->getNumberOfInputParameters(); j++)
            INFO("* IN  %d: %s", j, pcs->getInputParameterName(j).c_str());
        for (size_t j=0; j<pcs->getNumberOfOutputParameters(); j++)
            INFO("* OUT %d: %s", j, pcs->getOutputParameterName(j).c_str());
        ogs6fem->list_pcs.insert(pcs_name, pcs);
        auto opPCSList = opOgs.get_child_optional("processList");
        const ptree* opPCS = nullptr;
        if (opPCSList) {
            auto itr_range = opPCSList->equal_range("process");
            for (auto v=itr_range.first; v!=itr_range.second; ++v) {
                if (auto o = v->second.get_optional<std::string>("name")) {
                    if (o.get()==pcs_name) {
                        opPCS = &v->second;
                        break;
                    }
                }
            }
        }
        if (opPCS==nullptr) INFO("* Process option not found.");
        bool isPcsReady = pcs->initialize(opPCS!=nullptr ? *opPCS : opOgs);
        if (!isPcsReady) {
            ERR("***Error while setting up processes");
            return 0;
        }
    }

    INFO("->Setting time stepping...");
    TimeSteppingControllerWithOutput timestepping(&ogs6fem->outController);
    timestepping.setTransientSystem(*_cpl_system);

    //TODO the following calculation should be done in TimeSteppingController
    double t_start = std::numeric_limits<double>::max();
    double t_end = -1 * std::numeric_limits<double>::max();
    if (ogs6fem->list_tim.size() > 0) {
        for (size_t i=0; i<ogs6fem->list_tim.size(); i++) {
            t_start = std::min(t_start, ogs6fem->list_tim[i]->begin());
            t_end = std::max(t_end, ogs6fem->list_tim[i]->end());
        }
    } else {
        INFO("Time step configuration not found.");
        t_start = 0.0;
        t_end = 1.0;
    }

    INFO("->Outputting the initial values...");
    ogs6fem->outController.outputData(NumLib::TimeStep(t_start));

    //-------------------------------------------------------------------------
    // Run simulation
    //-------------------------------------------------------------------------
    INFO("->Simulation is ready! start=%f, end=%f", t_start, t_end);
    mpi.barrier();

    BaseLib::RunTime runTime;
    runTime.start();

    timestepping.setBeginning(t_start); //TODO really need this? start, end is already given in timestep function
    size_t n_timesteps = timestepping.solve(t_end);

    INFO("");
    INFO("->Simulation is finished.");
    INFO("#############################################################");
    INFO("*** Summary of this simulation");
    INFO("total time step : %d", n_timesteps);
    INFO("elapsed time   : %g sec", runTime.elapsed());
    INFO("#############################################################");
    INFO("");

    return 0;
}

} //end ogs6
