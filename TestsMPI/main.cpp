
#include <iostream>
#include <stdio.h>
#include "mpi.h"
#include "lis.h"

#include "MeshLib/Tools/MeshGenerator.h"

#include "DiscreteLib/LisDiscreteSystem.h"
#include "DiscreteLib/LisMPIDiscreteVector.h"
#include "DiscreteLib/ogs5/par_ddc_group.h"
#include "DiscreteLib/SparsityBuilder.h"

using namespace DiscreteLib;
using namespace MathLib;
using namespace MeshLib;

int nprocs = 0;
int my_rank = 0;

struct DiscreteExample1
{
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;

    DiscreteExample1()
    {
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    void setLocalDirichletBC(const Base::BidirectionalMap<size_t, size_t> &map_global2localNodeId, std::vector<size_t> &local_dirichlet_bc_id, std::vector<double> &local_dirichlet_bc_value)
    {
        for (size_t i=0; i<list_dirichlet_bc_id.size(); i++) {
            if (map_global2localNodeId.countInA(list_dirichlet_bc_id[i])>0) {
                size_t local_id = map_global2localNodeId.mapAtoB(list_dirichlet_bc_id[i]);
                local_dirichlet_bc_id.push_back(local_id);
                local_dirichlet_bc_value.push_back(list_dirichlet_bc_value[i]);
            }
        }
    }


    class TestElementAssembler : public IElemenetLocalAssembler
    {
        Matrix<double> _m;
    public:
        TestElementAssembler()
        {
            _m.resize(4,4);
            _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
            _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
            _m(2,2) = 4.0; _m(2,3) = -1.0;
            _m(3,3) = 4.0;
            for (size_t i=0; i<4; i++)
                for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);
            _m *= 1.e-11/6.0;
        }
        void assembly(MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs)
        {
            (*eqs.getA()) = _m;
        }
    };
};

void initialize(int argc, char *argv[])
{
    LisSolver::initialize(argc, argv);

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#else
    nprocs  = 1;
    my_rank = 0;
#endif
    if( my_rank==0 )
    {
        printf("\n");
        printf("number of processes = %d\n",nprocs);
    }
}

void finalize()
{
    LisSolver::finalize();
}

int testOGS5(int argc, char *argv[])
{
    OGS5::CPARDomainGroup *master;

    return 0;
}

void testLISSystem(int argc, char *argv[])
{
    // create test data
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    MeshLib::IMesh *org_msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    Base::BidirectionalMap<size_t, size_t> msh_node_id_mapping;
    std::set<size_t> ghost_nodes;
    MeshLib::IMesh* local_msh;
    if (my_rank==0) {
        int dom1_eles[] = {0, 1, 2, 3};
        int dom1_ghost_nodes[] = {5, 6, 7, 8};
        std::vector<size_t> dom1_e(dom1_eles, dom1_eles+4);
        MeshGenerator::generateSubMesh(*org_msh, dom1_e, local_msh, msh_node_id_mapping);
        for (size_t i=0; i<4; i++) {
            ghost_nodes.insert(msh_node_id_mapping.mapAtoB(dom1_ghost_nodes[i]));
        }
    } else {
        int dom2_eles[] = {1, 2, 3};
        int dom2_ghost_nodes[] = {0, 1, 2, 3, 4};
        std::vector<size_t> dom2_e(dom2_eles, dom2_eles+3);
        MeshGenerator::generateSubMesh(*org_msh, dom2_e, local_msh, msh_node_id_mapping);
        for (size_t i=0; i<5; i++) {
            ghost_nodes.insert(msh_node_id_mapping.mapAtoB(dom2_ghost_nodes[i]));
        }
    }
    MathLib::LisMPILinearEquation lis;
    lis.getOption().ls_method = LIS_option::CG;
    lis.getOption().ls_precond = LIS_option::NONE;


    // create linear system
    DofMapManager dof;
    LisMPILinearSystem<SparsityBuilderFromNodeConnectivity> linear_eq(*local_msh, lis, dof);
    linear_eq.setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
    linear_eq.construct(ElementBasedAssembler(ele_assembler));
    //linear_eq->getLinearEquation()->printout();
    linear_eq.solve();
}

void testLISVector(int argc, char *argv[])
{
    if (my_rank==0) 
        std::cout << "# testLISVector start" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    int n_lobal = 0;
    int n_global = 10;
    std::cout << my_rank << ": create vector" << std::endl;
    std::vector<int> list_ghost(2);
    if (my_rank==0) {
        list_ghost[0] = 5;
        list_ghost[1] = 6;
    } else {
        list_ghost[0] = 0;
        list_ghost[1] = 1;
    }
    LisMPIDiscreteVector v(MPI_COMM_WORLD, n_lobal, n_global, list_ghost);
    std::cout << my_rank << ": set vector" << std::endl;
    for (int i=v.getRangeBegin(); i<v.getRangeEnd(); i++) {
        v.global(i) = i;
    }
    std::cout << my_rank << ": finish vector" << std::endl;
    v.finishUpdate();

    std::cout << my_rank << ": output" << std::endl;
    v.print();
    std::vector<double> x(v.getGlobalSize());
    v.getGlobal(&x[0]);
    if (my_rank==0) {
        std::cout << "x=" << std::endl;
        for (size_t i=0; i<x.size(); i++)
            std::cout << x[i] << " ";
        std::cout << std::endl;
    }
    std::cout << my_rank << ": ghost elements (" << v.getNumberOfGhostElements() << ")" << std::endl;
    for (size_t i=0; i<v.getNumberOfGhostElements(); i++) {
        std::cout << my_rank << ": [" << v.getGhostElementId(i) << "]=" << v.global(v.getGhostElementId(i)) << std::endl;
    }
    std::cout << my_rank << ": end" << std::endl;
}

int testLIS(int argc, char *argv[])
{
    // input parameters
    if( argc < 3 ) {
        if( my_rank==0 ) printf("Usage: test5 n gamma [options]\n");
        return -1;
    }

    int gn  = atoi(argv[1]);
    double gamma  = atof(argv[2]);

    // solve
    LisSolver lis;
#if 0
    SparseTableCRS<int>* crs = lis.createCRS(0, gn);
    int is = 0;
    int ie = 0;
    lis.getRange(is,ie);

    int k = 0;
    crs->row_ptr[0] = 0;
    for (int ii=is;ii<ie;ii++) {
        int jj=0;
        if( ii>1 )    { jj = ii - 2; crs->col_idx[k] = jj; crs->data[k++] = gamma;}
        if( ii<gn-1 ) { jj = ii + 1; crs->col_idx[k] = jj; crs->data[k++] = 1.0;}
        crs->col_idx[k] = ii; crs->data[k++] = 2.0;
        crs->row_ptr[ii-is+1] = k;
    }
#else
    lis.createDynamic(0, gn);
    int is = 0;
    int ie = 0;
    lis.getRange(is,ie);
    std::cout << my_rank << ": is=" << is << ", ie=" << ie << std::endl;
    for (int i=is;i<ie;i++) {
        if( i>1   )  lis.addA(i,i-2,gamma);
        if( i<gn-1 ) lis.addA(i,i+1,1.0);
        lis.addA(i,i,2.0);
    }
#endif

    lis.assembleMatrix();
    size_t iu = lis.createVector();
    lis.setVectorAll(iu, 1.0);
    lis.matvecToRHS(iu);
    lis.destroyVector(iu);
    lis.solve();

    return 0;
}

int main(int argc, char *argv[])
{
    initialize(argc, argv);
    testLISVector(argc, argv);
    //testLIS(argc, argv);
    finalize();
#if 0
    MPI_Init(&argc, &argv);

    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank==0) {
        std::cout << "MPI program started with " << size << " processes " << std::endl;
        std::cout << std::endl;
        std::cout << "Who are you?" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Hello, I am process id " << rank  << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0) {
        std::cout << std::endl << "Talk each other!" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int say[] = {1, 2, 3, 4, 5, 6, 7};
    int hear[10] = {};
    if (rank%2==0) {
        int talk_to = (rank+1)%size;
        int hear_from = (rank+1)%size;
        int n_say = rank;
        int n_hear = rank+1;
        MPI_Send(say, n_say, MPI_INT, talk_to, 0, MPI_COMM_WORLD);
        std::cout << rank <<  " says hello to " << talk_to << std::endl;
        MPI_Status status;
        MPI_Recv(hear, n_hear, MPI_INT, talk_to, 0, MPI_COMM_WORLD, &status);
        std::cout << rank << " hear hello from " << talk_to << std::endl;
    } else {
        int talk_to = (rank-1)%size;
        int hear_from = (rank-1)%size;
        int n_say = rank;
        int n_hear = rank-1;
        MPI_Status status;
        MPI_Recv(hear, n_say, MPI_INT, talk_to, 0, MPI_COMM_WORLD, &status);
        std::cout << rank << " hear from " << talk_to << std::endl;
        MPI_Send(say, n_hear, MPI_INT, talk_to, 0, MPI_COMM_WORLD);
        std::cout << rank <<  " says hello to " << talk_to << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) {
        std::cout << std::endl << "I distribute this to you!" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int dist_v = 123;
    MPI_Bcast(&dist_v, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::cout << rank <<  " received a number " << dist_v << std::endl;

    MPI_Finalize();
#endif
    return 0;
}

