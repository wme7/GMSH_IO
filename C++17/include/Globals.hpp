#ifndef GLOBALS
#define GLOBALS

#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <map>

// local libs
#include "MdimArray.hpp"
#include "NpyArray.hpp"

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

        // Control Objects
        MArray<size_t,2> PointsToV;
        MArray<size_t,2> LinesToV;
        MArray<size_t,2> FacesToV;

        // Member functions
        void print();
        void printLocal();
        // get mesh info
        size_t get_D() const {return spatialDimension;};
        size_t get_Nfaces() const {return Nfaces;};
        size_t get_nPartitions() const {return num_partitions;};
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

#endif