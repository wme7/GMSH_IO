#ifndef GMSH_PARSER_V2
#define GMSH_PARSER_V2

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

// Get a sub-string from the main string-buffer using two unique delimiters:
std::string extractBetween(
    const std::string &buffer,
    const std::string &start_delimiter,
    const std::string &stop_delimiter)
{
    unsigned first_delim_pos = buffer.find(start_delimiter);
    unsigned end_pos_of_first_delim = first_delim_pos + start_delimiter.length();
    unsigned last_delim_pos = buffer.find(stop_delimiter);

    return buffer.substr(end_pos_of_first_delim,
        last_delim_pos - end_pos_of_first_delim);
}

// Build convention map of Boundary conditions for ParadigmS:
std::map<std::string,int> get_BC_type()
{
    std::map<std::string,int> BC_type;

    // Initialize map of BCs
    BC_type["BCfile"]=0;
    BC_type["free"]=1;
    BC_type["wall"]=2;
    BC_type["outflow"]=3;
    BC_type["imposedPressure"]=4;
    BC_type["imposedVelocities"]=5;
    BC_type["axisymmetric_y"]=6;
    BC_type["axisymmetric_x"]=7;
    BC_type["BC_rec"]=10;
    BC_type["free_rec"]=11;
    BC_type["wall_rec"]=12;
    BC_type["outflow_rec"]=13;
    BC_type["imposedPressure_rec"]=14;
    BC_type["imposedVelocities_rec"]=15;
    BC_type["axisymmetric_y_rec"]=16;
    BC_type["axisymmetric_x_rec"]=17;
    BC_type["piston_pressure"]=18;
    BC_type["piston_velocity"]=19;
    BC_type["recordingObject"]=20;
    BC_type["recObj"]=20;
    BC_type["piston_stress"]=21;

    return BC_type;
}

// Element structure to be acquired:
struct element
{
    std::vector<size_t> EToV;
    std::vector<int> phys_tag;
    std::vector<int> geom_tag;
    std::vector<int> part_tag;
    std::vector<int> Etype;
};

// Split and input strings and return a vector with all double values
std::vector<double> str2double(const std::string &str)
{
    std::stringstream ss(str);
    std::istream_iterator<double> begin(ss);
    std::istream_iterator<double> end;
    std::vector<double> values(begin,end);
    return values;
}

// Split and input strings and return a vector with all size_t values
std::vector<size_t> str2size_t(const std::string &str)
{
    std::stringstream ss(str);
    std::istream_iterator<size_t> begin(ss);
    std::istream_iterator<size_t> end;
    std::vector<size_t> values(begin,end);
    return values;
}

// Get nodes from block system
MArray<double,2> get_nodes(const std::string &Nodes, const size_t &numNodes) 
{
    // Allocate output
    MArray<double,2> V({numNodes,3},0); // [x(:),y(:),z(:)]

    // Node counter
    size_t nodeTag;

    std::stringstream buffer(Nodes);
    std::string line; 
    std::getline(buffer, line); // l = 1;
    // Read nodes blocks:  
    for (size_t i=0; i<numNodes; i++)
    {
        std::getline(buffer, line);
        std::stringstream stream(line);
        stream >> nodeTag >> V(i,0) >> V(i,1) >> V(i,2);
    }
    return V;
}

int GMSHparserV2(std::string mesh_file) 
{
    //-------------------------------------
    // - Global Parameters of the parser -
    //-------------------------------------

    // Nominal dimension of the meshed geometry
    size_t phys_DIM = 0; // initial guess

    // Is the mesh a single partition?
    bool single_domain = true; // Initial guess

    // Index correction factor
    size_t one = 0; // if we set one=0, one recovers the original indexes of GMSH.

    // Print comments 
    bool DEBUG = false; // print all all warnings and stage messages

    // Build Maps variables
    std::map<size_t, std::string> phys2names;
    std::map<size_t, int> point2phys, point2part, point2geom;
    std::map<size_t, int> curve2phys, curve2part, curve2geom;
    std::map<size_t, int> surf2phys, surf2part, surf2geom;
    std::map<size_t, int> volm2phys, volm2part, volm2geom;

    // Last by not least, open meshfile:    
    std::ifstream file(mesh_file);

    //-------------------------------------
    // -- READ GMSH V2.2 FORMAT --
    //-------------------------------------   
    if (file)
    {   
        // Fill buffer
        std::stringstream buffer;
        buffer << file.rdbuf();
        file.close();

        // Initialized sub-buffers
        std::string MeshFormat    = extractBetween(buffer.str(),"$MeshFormat\n","\n$EndMeshFormat");
        std::string PhysicalNames = extractBetween(buffer.str(),"$PhysicalNames\n","\n$EndPhysicalNames");
        std::string Entities      = extractBetween(buffer.str(),"$Entities\n","\n$EndEntities");
        std::string PartEntities  = extractBetween(buffer.str(),"$PartitionedEntities\n","\n$EndPartitionedEntities");
        std::string Nodes         = extractBetween(buffer.str(),"$Nodes\n","\n$EndNodes");
        std::string Elements      = extractBetween(buffer.str(),"$Elements\n","\n$EndElements");

        // Sanity check
        if (  MeshFormat.empty() ) {std::cout << " Error - Wrong File Format!" << std::endl; std::exit(-1);}
        if (PhysicalNames.empty()) {std::cout << " Error - No Physical names!" << std::endl; std::exit(-1);}
        if (   Entities.empty()  ) {std::cout << " Error - No Entities found!" << std::endl; std::exit(-1);}
        if (    Nodes.empty()    ) {std::cout << " Error - Nodes are missing!" << std::endl; std::exit(-1);}
        if (   Elements.empty()  ) {std::cout << " Error - No Elements found!" << std::endl; std::exit(-1);}

        // Is it a single or partitioned domain?
        if ( PartEntities.empty()) {single_domain = true;} else {single_domain = false;}
        if(DEBUG) std::cout << "Partitioned domain detected" << std::endl;

        // For reading the buffer line by line
        std::string line;

        /**************************/
        // 1. Read MeshFormat
        /**************************/
        std::stringstream buffer_MF;    
        buffer_MF << MeshFormat;

        double version = 1.0; size_t format=0, size=0;
        buffer_MF >> version >> format >> size;
        std::cout << "Mesh version " << version << ", Binary " << format << ", endian " << size << std::endl;

        // Sanity check
        if (version != 2.2)
        {
            std::cout << " Error - Expected mesh format v4.1" << std::endl; 
            std::exit(-1);
        }
        if (format != 0)
        {
            std::cout << " Error - Binary file not allowed" << std::endl; 
            std::exit(-1);
        }
        // Clear Buffer
        buffer_MF.str(std::string());  // this is equivalent to: buffer.str("");

        /**************************/
        // 2. Read PhysicalNames
        /**************************/
        std::stringstream buffer_PN;
        buffer_PN << PhysicalNames;

        size_t num_physical_groups = 0;
        buffer_PN >> num_physical_groups;

        if(DEBUG) std::cout << "num_physical_groups: " << num_physical_groups << std::endl;

        for (size_t i=0; i<num_physical_groups; ++i)
        {
            // Read an entire line of the PhysicalNames section.
            std::getline(buffer_PN, line);
            size_t phys_dim, phys_id;
            std::string phys_name;
            buffer_PN >> phys_dim >> phys_id >> phys_name;

            // get rid of the quotes characters from the phys_name
            phys_name.erase(std::remove(phys_name.begin(), phys_name.end(), '"'), phys_name.end());

            // create a map between phys_id to phys_name.
            phys2names[phys_id] = phys_name;

            // Search for the maximun dimension of entities in the mesh
            phys_DIM = phys_DIM > phys_dim ? phys_DIM : phys_dim;

            if(DEBUG) std::cout << " Physical Name: " << phys_name << std::endl;
            if(DEBUG) std::cout << " Physical ID: " << phys_id << std::endl;
            if(DEBUG) std::cout << " Entity Dim: " << phys_dim << std::endl;
        }
        // Clear Buffer
        buffer_PN.str(std::string());

        /**************************/
        // 3. Read Nodes
        /**************************/
        std::stringstream buffer_N;
        buffer_N << Nodes;

        size_t numNodes = 0;
        buffer_N >> numNodes;

        if(DEBUG) std::cout << " numNodes: " << numNodes << std::endl;

        //V = get_nodes(buffer_N,runNodesBlocks,numNodes);
        MArray<double,2> V = get_nodes(Nodes,numNodes);

        // Report the number for nodes found
        std::cout << "Total vertices found = " << numNodes << std::endl;

        // Clear Buffer
        buffer_N.str(std::string());

        /**************************/
        // 4. Read Elements
        /**************************/
        std::stringstream buffer_E;
        buffer_E << Elements;

        size_t numElements  = 0;
        buffer_E >> numElements;

        if(DEBUG) std::cout << " numElements: " << numElements << std::endl;
        std::getline(buffer_E, line); // if ok, go to next line ..

        // Allocate space for Elements and their data
        element PE; // point elements
        element LE; // line elements
        element SE; // surface elements
        element VE; // volume elements

        // Element counters
        size_t e1 = 0; // Lines Element counter
        size_t e2 = 0; // Triangle Element counter
        size_t e4 = 0; // Tetrahedron Element counter
        size_t e15= 0; // Point Element counter

        // Read elements serially
        for (size_t i=0; i<numElements; i++) 
        {
            // Read Block parameters
            std::getline(buffer_E, line);
        }
        std::cout << "Total point-elements found = " << e15 << std::endl;
        std::cout << "Total curve-elements found = " << e1 << std::endl;
        std::cout << "Total surface-elements found = " << e2 << std::endl;   
        std::cout << "Total volume-elements found = " << e4 << std::endl;
        // Sanity check
        if (numElements != (e15+e1+e2+e4)) {
            std::cout << "Total number of elements missmatch!"<< std::endl;
            std::exit(-1);
        }
        // Clear Buffer
        buffer_E.str(std::string());

        /**************************/
        // 5. Save parsed arrays
        /**************************/

        // Save to Numpy Array
        //cnpy::npy_save("FEMmesh.npy",&PE.EToV[0],{e15,1},"w");
        //cnpy::npy_save("FEMmesh.npy",&LE.EToV[0],{ e1,2},"a");
        //cnpy::npy_save("FEMmesh.npy",&SE.EToV[0],{ e2,3},"a");
        //cnpy::npy_save("FEMmesh.npy",&VE.EToV[0],{ e4,4},"a");

        /**************************/
        // 6. Print parsed arrays
        /**************************/

        // Save to MdimArrays
        // MArray<size_t,2> PEToV({e15,1},PE.EToV); PEToV.print();
        // MArray<size_t,2> LEToV({ e1,2},LE.EToV); LEToV.print();
        // MArray<size_t,2> SEToV({ e2,3},SE.EToV); SEToV.print();
        // MArray<size_t,2> VEToV({ e4,4},VE.EToV); VEToV.print();

    } else {
        std::cout << "ERROR: Could not open file: " << mesh_file << std::endl; 
        std::exit(-1);
    }
    // If everything goes well ...
    return 0;
}

#endif