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
    // Test the parser V2 and V4 (similarly to Matlab prototype)
    GMSHparserV2("../../meshes/rectangle_v2.msh");
    GMSHparserV2("../../meshes/cuboid_v2.msh");
    GMSHparserV4("../../meshes/rectangle_v4.msh");
    GMSHparserV4("../../meshes/cuboid_v4.msh");

    // Note: The verification and visualization of the mesh is done with a script in python3

    // If everything is normal, then ..
    return 0;
}