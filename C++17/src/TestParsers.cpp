//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//  Test GMSH parsers for msh-files written in format version v2.2 and V4.1
//
//      Coded by Manuel A. Diaz @ Pprime | Univ-Poitiers, 2022.01.21
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "Globals.hpp"
#include "GMSHparserV2.hpp"
#include "GMSHparserV4.hpp"

//int main(int argc, char *argv[]) 
int main() 
{
    // Load Objects
    //FEMmesh mesh;

    // if (argc<2)
    // { 
    //     std::cout << "WARNING" << std::endl;
    //     std::cout << "The msh file is missing" << std::endl;
    // } 
    // else
    // {
    //     std::string filename = argv[1];
    //     std::cout << "The msh file is :" << filename << std::endl;
        
    //     // Get mesh data and parameters
    //     //mesh = GMSHparser(parameters.mesh_name());
    //     GMSHparserV4(filename);
    // }
    // verify computation
    //mesh.print();

    GMSHparserV2("../../meshes/rectangle_v2.msh");
    //GMSHparserV2("../../meshes/cuboid_v2.msh");
    GMSHparserV4("../../meshes/rectangle_v4.msh");
    //GMSHparserV4("../../meshes/cuboid_v4.msh");

    // If everything is normal
    return 0;
}