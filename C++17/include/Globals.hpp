#ifndef GLOBALS
#define GLOBALS

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <algorithm>
#include "MdimArray.hpp"
#include "NpyArray.hpp"
#include <mpi.h>
//#include <omp.h>

// Objects required to start DG2D.hpp or DG3D.hpp

// --Parameters for simulation-- object
class Parameters
{
    //access
    public:

        // Basic Parameters
        size_t m_N;         // Polynomial base degree
        size_t m_Npara;     // Number of parameters
        size_t m_sv;        // spacing of the saving points
        size_t m_Nt_user;   // requested Number of (fixed) time steps
        size_t m_nb_sources;// Number of point (Gaussian) sources
        size_t m_nb_boundaries_sources;// Number of point BC (Gaussian) sources
        double m_CFL;       // Courant–Friedrichs–Lewy condition parameter ( CFL < 1.0 )
        double m_dt_user;   // requested (fixed) time step by the user.
        
        // Simulation, Solver and Mesh names
        std::string m_name;
        std::string m_solver;
        std::string m_mesher;
        std::string m_binary;
        std::string m_mesh_name;
        std::string m_SSfile;
        std::string m_BCSSfile;

        // IC, BCs, STs, PhysParams and output modes
        std::string m_initial_mode;
        std::string m_boundaries_mode;
        std::string m_sources_mode;
        std::string m_physical_parameter_mode;
        std::string m_xmf_mode;
        std::string m_ipynb_mode;
        std::string m_terminal_mode;
        std::string m_query_fields_mode;

        // Query max and min fields and their initial values
        int m_qmin_eq, m_qmax_eq, m_qsum_eq, m_qsum2_eq, m_qsum3_eq, m_qsum4_eq;
        double m_qmin_0, m_qmax_0, m_qsum_0, m_qsum2_0, m_qsum3_0, m_qsum4_0;
        int m_qIntensity_eq, m_qVelocity_eq, m_qSkewness_eq, m_qKurtosis_eq, m_qRMS_eq;
        double m_qIntensity_0, m_qVelocity_0, m_qSkewness_0, m_qKurtosis_0, m_qRMS_0;
        double m_tq_start, m_tq_end;

        // Physical Parameters of the medium of propagation
        std::vector<std::string> m_physical_parameter_names;
        std::vector<double> m_physical_parameter_values;

        // For Gaussian Initial Conditions:
        double m_Amp;
        double m_x0;
        double m_y0;
        double m_z0;
        double m_Sgm;
        double m_R;
        double m_theta;

        // For Imposed Boundary Conditions:
        double m_p0_BC;
        double m_u0_BC;
        double m_v0_BC;
        double m_w0_BC;
        double m_f0_BC;

        // Source parameters
        std::vector<double> m_A;
        std::vector<double> m_f0;
        std::vector<double> m_sigma;
        std::vector<double> m_xs;
        std::vector<double> m_ys;
        std::vector<double> m_zs;
        std::vector<double> m_T0;
        std::vector<double> m_phi;
        std::vector<double> m_sources_eq;

        // Simulation time
        double m_t_start;
        double m_t_end;

        // recording objects parameters
        std::string m_rec_mode;
        int m_rec_eq;

        // Member functions
        void print();
        void printSources();
        size_t polyOrder() const {return m_N;};
        std::string name() const {return m_name;};
        std::string solver() const {return m_solver;};
        std::string mesh_name() const {return m_mesh_name;};  
        const char* get_TestName() const {return m_name.c_str();};
        const char* get_SolverName() const {return m_solver.c_str();};
        const char* get_MeshfileName() const {return m_mesh_name.c_str();};
};

//=============================================================
void Parameters::print()
{
    std::cout << "Parametres de la simulation " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "      Name  : " << m_name <<std::endl;
    std::cout << "      Mesh  : " << m_mesh_name <<std::endl;
    std::cout << "    Solver  : " << m_solver <<std::endl;
    std::cout << "     P-deg  : " << m_N <<std::endl;
    std::cout << "   IC mode  : " << m_initial_mode <<std::endl;
    std::cout << "   BC mode  : " << m_boundaries_mode <<std::endl;
    std::cout << "   SS mode  : " << m_sources_mode <<std::endl;
    std::cout << "   PP mode  : " << m_physical_parameter_mode <<std::endl;
    std::cout << "  nSS input : " << m_nb_sources <<std::endl;
    std::cout << "  nPP input : " << m_physical_parameter_names.size() <<std::endl;
    std::cout << "   t_start  : " << m_t_start <<std::endl;
    std::cout << "   t_output : " << m_t_end <<std::endl;
    std::cout << "===========================================" << std::endl;
}

void Parameters::printSources()
{
    if (m_sources_mode != "false")
    {
        std::cout << "Sources in the simulation " << std::endl;
        std::cout << "===========================================" << std::endl;
        std::cout << " Number of sources :\t" << m_nb_sources <<std::endl;
        std::cout << "   A :\t"; for (auto m:m_A) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << "  f0 :\t"; for (auto m:m_f0) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << " sgm :\t"; for (auto m:m_sigma) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << "  xc :\t"; for (auto m:m_xs) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << "  ys :\t"; for (auto m:m_ys) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << "  zs :\t"; for (auto m:m_zs) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << "  T0 :\t"; for (auto m:m_T0) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << "  Eq :\t"; for (auto m:m_sources_eq) std::cout <<m<< "\t"; std::cout << std::endl;
        std::cout << "===========================================" << std::endl;
    }
}

// --FEM Mesh (or Skeleton Mesh) data-- Object
class FEMmesh
{
    //access
    public:
        // Fundamental parameteres of the discrete domain
        size_t spatialDimension=0; // Domain Spatial Dimension
        size_t iP=0; // rank/partition ID number (assumed to be ROOT partition as initial guess)
        size_t num_partitions=0; // Number of partitions in the mesh
        std::vector<size_t> KP; // map of domain-elements per partition
        std::vector<size_t> KBP; // map of boundary-elements per partition
        std::vector<size_t> CO_Points_P; // map of control points per partition
        std::vector<size_t> CO_Lines_P; // map of control lines per partition
        std::vector<size_t> CO_Faces_P; // map of control faces per partition

        // Global mesh parameters
        size_t K=0; // number of domain elements
        size_t KB=0; // number of boundary elements
        size_t Nfaces=0; // number faces of domain elements
        size_t num_nodes=0; // number of nodes in the domain

        // Global control objects inside the mesh
        size_t CO_Points=0;
        size_t CO_Lines=0;
        size_t CO_Faces=0;

        // Global Mesh information
        MArray<double,2> V;
        MArray<double,2> EToVc;
        MArray<size_t,2> EToV;
        MArray< int, 2 > BEToV;
        MArray<size_t,2> EToP;
        MArray<size_t,2> BEToP;
        std::vector<int> Etype;
        MArray< int, 2 > EToE;
        MArray<size_t,2> EToF;
        MArray< int, 2 > EToBE;

        // Control Objects
        MArray<size_t,2> PointsToV;
        MArray<size_t,2> LinesToV;
        MArray<size_t,2> FacesToV;

        // Local Mesh parameters
        size_t Ki=0; // number of elements in the local domain (or mesh)
        size_t KBi=0; // number of boundary elements in the local mesh
        size_t KBufferi=0; // number of buffer elements in the local mesh

        // Local control objects inside the mesh
        size_t CO_Points_i=0;
        size_t CO_Lines_i=0;
        size_t CO_Faces_i=0;

        // Local to Global elements indexes
        std::vector<int> LToG;

        // Local Mesh information
        MArray<double,2> L_EToVc;
        MArray<size_t,2> L_EToV;
        MArray< int, 2 > L_BEToV;
        MArray< int, 2 > L_BufferEToV;
        MArray< int, 2 > L_CommPattern;
        std::vector<int> L_Etype;
        MArray< int, 2 > L_EToE;
        MArray<size_t,2> L_EToF;
        MArray< int, 2 > L_EToBE;

        // Local Control Objects
        MArray<size_t,2> L_PointsToV;
        MArray<size_t,2> L_LinesToV;
        MArray<size_t,2> L_FacesToV;
        
        std::vector<std::vector<int>> L_PointsToE;
        std::vector<std::vector<int>> L_LinesToE;
        std::vector<std::vector<int>> L_FacesToE;
        std::vector<std::vector<int>> L_PointsToNode;
        std::vector<std::vector<int>> L_LinesToEdge;
        std::vector<std::vector<int>> L_FacesToFace;

        // time step scale
        double dt_scale = 0; // Global dt_scale
        double L_dt_scale = 0; // Local dt_scale

        // rmin = abs(rLGL(0)-rLGL(1));
        double rmin = 0; // Global rmin
        double L_rmin = 0; // Local rmin

        // Member functions
        void print();
        void printLocal();
        // get mesh info
        size_t get_D() const {return spatialDimension;};
        size_t get_Nfaces() const {return Nfaces;};
        size_t get_nPartitions() const {return num_partitions;};
        // get local info
        size_t get_K() const {return Ki;};
        size_t get_KB() const {return KBi;};
        size_t get_KBuffer() const {return KBufferi;};
        size_t get_CO_Points() const {return CO_Points_i;};
        size_t get_CO_Lines() const {return CO_Lines_i;};
        size_t get_CO_Faces() const {return CO_Faces_i;};
        // get global info
        size_t get_global_K() const {return K;};
        size_t get_global_KB() const {return KB;};
        size_t get_global_nNodes() const {return num_nodes;};
        size_t get_global_CO_Points() const {return CO_Points;};
        size_t get_global_CO_Lines() const {return CO_Lines;};
        size_t get_global_CO_Faces() const {return CO_Faces;};
};

//=============================================================
void FEMmesh::print()
{
    std::cout << "\n Global parameters of mesh stored in rank " << iP <<std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << "  D :" << spatialDimension << std::endl;
    std::cout << " nV :" << num_nodes << std::endl;
    std::cout << " nP :" << num_partitions << std::endl;
    std::cout << "  K :" << K << std::endl;
    std::cout << " KB :" << KB << std::endl;
    std::cout << " KP :"; for (auto k:KP) std::cout << k << " "; std::cout << std::endl;
    std::cout << " KBP:"; for (auto k:KBP) std::cout << k << " "; std::cout << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << "\n Control objects in mesh stored in rank " << iP <<std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << " CO_Points :" << CO_Points << std::endl;
    std::cout << " CO_Lines  :" << CO_Lines << std::endl;
    std::cout << " CO_Faces  :" << CO_Faces << std::endl;
    std::cout << " CO_Points_P :"; for (auto co:CO_Points_P) std::cout <<co<< " "; std::cout << std::endl;
    std::cout << " CO_Lines_P  :"; for (auto co:CO_Lines_P) std::cout <<co<< " "; std::cout << std::endl;
    std::cout << " CO_Faces_P  :"; for (auto co:CO_Faces_P) std::cout <<co<< " "; std::cout << std::endl;
    std::cout << "===================================================" << std::endl;
    if (V.total_size()>0)       cnpy::npz_save("FEM_Geometry.npz",  "V"  ,  V.data()  ,  V.dims()  ,"w");
    if (EToV.total_size()>0)    cnpy::npz_save("FEM_Geometry.npz","EToV" ,EToV.data() ,EToV.dims() ,"a");
    if (BEToV.total_size()>0)   cnpy::npz_save("FEM_Geometry.npz","BEToV",BEToV.data(),BEToV.dims(),"a");
    if (EToP.total_size()>0)    cnpy::npz_save("FEM_Geometry.npz","EToP" ,EToP.data() ,EToP.dims() ,"a");
    if (BEToP.total_size()>0)   cnpy::npz_save("FEM_Geometry.npz","BEToP",BEToP.data(),BEToP.dims(),"a");
    if (PointsToV.total_size()>0) cnpy::npz_save("FEM_Geometry.npz","PointsToV",PointsToV.data(),PointsToV.dims(),"a");
    if (LinesToV.total_size()>0)  cnpy::npz_save("FEM_Geometry.npz","LinesToV" ,LinesToV.data() ,LinesToV.dims() ,"a");
    if (FacesToV.total_size()>0)  cnpy::npz_save("FEM_Geometry.npz","FacesToV" ,FacesToV.data() ,FacesToV.dims() ,"a");
}

//======================FEM_=======================================
void FEMmesh::printLocal()
{
    std::cout << "\n Local parameters of partitioned mesh stored in rank " << iP <<std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "  D :" << spatialDimension << std::endl;
    std::cout << "  K :" << Ki << std::endl;
    std::cout << " KB :" << KBi << std::endl;
    std::cout << " KBuffer :" << KBufferi << std::endl;
    std::cout << " dt_scale :" << L_dt_scale << std::endl;
    std::cout << " rmin :" << L_rmin << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "\n Control objects in mesh stored in rank " << iP <<std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << " CO_Points :" << CO_Points_i << std::endl;
    std::cout << " CO_Lines :" << CO_Lines_i << std::endl;
    std::cout << " CO_Faces :" << CO_Faces_i << std::endl;
    std::cout << "==================================================" << std::endl;
    if (V.total_size()>0)       cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz",   "V"   ,   V.data()   ,   V.dims()    ,"w");
    if (LToG.size()>0)          cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz",  "LToG" ,  LToG.data() , {LToG.size()} ,"a");
    if (L_EToV.total_size()>0)  cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_EToV" ,L_EToV.data() , L_EToV.dims() ,"a");
    if (L_BEToV.total_size()>0) cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_BEToV",L_BEToV.data(), L_BEToV.dims(),"a");
    if (L_EToE.total_size()>0)  cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_EToE" ,L_EToE.data() , L_EToE.dims() ,"a");
    if (L_EToF.total_size()>0)  cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_EToF" ,L_EToF.data() , L_EToF.dims() ,"a");
    if (L_Etype.size()>0)       cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_Etype",L_Etype.data(),{L_Etype.size()},"a");
    if (L_PointsToV.total_size()>0)   cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_PointsToV" ,L_PointsToV.data() ,L_PointsToV.dims(),"a");
    if (L_LinesToV.total_size()>0)    cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_LinesToV"  ,L_LinesToV.data()  ,L_LinesToV.dims() ,"a");
    if (L_FacesToV.total_size()>0)    cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_FacesToV"  ,L_FacesToV.data()  ,L_FacesToV.dims() ,"a");
    if (L_BufferEToV.total_size()>0)  cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_BufferEToV" ,L_BufferEToV.data() ,L_BufferEToV.dims() ,"a");
    if (L_CommPattern.total_size()>0) cnpy::npz_save("FEM_geometry_d"+std::to_string(iP)+".npz","L_CommPattern",L_CommPattern.data(),L_CommPattern.dims(),"a");
}

//=============================================================

#endif