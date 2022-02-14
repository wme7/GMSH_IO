#ifndef GMSH_PARSER_V4
#define GMSH_PARSER_V4

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

// Get entity Tag and its associated physical Tag.
std::tuple<size_t, int> get_entity(const std::string &line, const size_t &Idx)
{
    size_t entityTag;
    size_t numPhysicalTags; 
    int physicalTag;

    // get line data
    auto vector = str2double(line);

    switch (Idx) {
        case 1: // case for nodes
            // 1. get entityTag
            entityTag = int(vector[0]);

            // 2. get entity coordiantes // not needed
            // ignore indexes 1, 2, 3

            // 3. get physical tag associated
            numPhysicalTags = int(vector[4]);
            if (numPhysicalTags==0) {
                physicalTag = -1; // set a negative tag!
            } else {
                physicalTag = int(vector[5]);
            }
            break;
        default: // otherwise
            // 1. get entityTag
            entityTag = int(vector[0]);

            // 2. get entity coordiantes // not needed
            // ignore indexes 1, 2, 3, 4, 5, 6

            // 3. get physical tag associated
            numPhysicalTags = int(vector[7]);
            if (numPhysicalTags==0) {
                physicalTag = -1; // set a negative tag!
            } else {
                physicalTag = int(vector[8]);
            }
            // 4. get tags of subentities. // not needed
            break;
    }
    return {entityTag, physicalTag};
}

// Get partitioned entity Tags and its associated physical Tag.
std::tuple<size_t, size_t, int, int> get_partitionedEntity(const std::string &line, const size_t &Idx)
{
    size_t entityTag, parentTag; 
    size_t numPhysicalTags, numPartitionsTags; 
    int partitionTag, physicalTag;

    // get line data
    auto vector = str2double(line);

    switch (Idx) {
        case 1: // case for nodes
            // 1. get entityTag
            entityTag = int(vector[0]);

            // 2. get parent dimension and tag
            //parentDim = int(vector[1]); // not needed
            parentTag = int(vector[2]);

            // 3. get partitions tags
            numPartitionsTags = int(vector[3]);
            if (numPartitionsTags > 1) { // --> mark it as an interface element!
                partitionTag = -1;
            } else {
                partitionTag = int(vector[4]);
            }

            // 4. get entity coordiantes // not needed
            // ignore indexes 5, 6, 7

            // 5. get physical tag associated
            numPhysicalTags = int(vector[7+numPartitionsTags]);
            if (numPhysicalTags==0) {
                physicalTag = -1; // set a negative tag!
            } else {
                physicalTag = int(vector[8+numPartitionsTags]);
            }
            break;
        default: //otherwise
            // 1. get entityTag
            entityTag = int(vector[0]);

            // 2. get parent dimension and tag
            //parentDim = int(vector[1]); // not needed
            parentTag = int(vector[2]);

            // 3. get partitions tags
            numPartitionsTags = int(vector[3]);
            if (numPartitionsTags > 1) { // --> mark it as an interface element!
                partitionTag = -1; // set a negative tag!
            } else {
                partitionTag = int(vector[4]);
            }

            // 4. get entity coordiantes // not needed
            // ignore indexes 5, 6, 7, 8, 9, 10

            // 5. get physical tag associated
            numPhysicalTags = int(vector[10+numPartitionsTags]);
            if (numPhysicalTags==0) {
                physicalTag = -1;
            } else {
                physicalTag = int(vector[12+numPartitionsTags]);
            }
            break;
    }
    return {entityTag, parentTag, partitionTag, physicalTag};
}

// Get nodes from block system
MArray<double,2> get_nodes(const std::string &Nodes, const size_t &numNodesBlocks, const size_t &numNodes) 
{
    // Allocate output
    MArray<double,2> V({numNodes,3}); V.set_value(0.0); // [x(:),y(:),z(:)]

    // Node counter
    size_t n=0;

    std::stringstream buffer(Nodes);
    std::string line; 
    std::getline(buffer, line); // l = 1;
    // Read nodes blocks:  (this can be done in parallel!)
    for (size_t i=0; i<numNodesBlocks; i++)
    {
        // update line counter, l = l+1;
        std::getline(buffer, line);
        std::stringstream hearder(line);

        // Read Block parameters
        size_t entityDim;  // not needed
        size_t entityTag;  // not needed
        size_t parametric; // not needed
        size_t numNodesInBlock;
        hearder >> entityDim >> entityTag >> parametric >> numNodesInBlock;

        // Read Nodes IDs
        size_t *nodeTag = new size_t[numNodesInBlock]; // nodeTag
        for (size_t i=0; i<numNodesInBlock; i++)
        {
            std::getline(buffer, line);
            std::stringstream stream(line);
            stream >> nodeTag[i];
        }
        // Read Nodes Coordinates
        for (size_t i=0; i<numNodesInBlock; i++)
        {
            std::getline(buffer, line);
            std::stringstream stream(line);
            stream >> V(n,0) >> V(n,1) >> V(n,2);
            n = n+1; // Update node counter
        }
        // Delete temporary new-arrays
        delete [] nodeTag;
    }
    return V;
}


int GMSHparserV4(std::string mesh_file) 
{
    //-------------------------------------
    // 0 - Global Parameters of the parser
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
    //  1 - READ GMSH V4.1 FORMAT
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

            // Sanity check
            if (version != 4.1)
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
        // Read PhysicalNames
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

        if (single_domain) {
            /**************************/
            // Read Entities
            /**************************/
            std::stringstream buffer_Ent;
            buffer_Ent << Entities;

            size_t numPoints  = 0;
            size_t numCurves  = 0;
            size_t numSurfaces= 0;
            size_t numVolumes = 0;
            buffer_Ent >> numPoints >> numCurves >> numSurfaces >> numVolumes;

            if(DEBUG) std::cout << " numPoints: "   << numPoints   << std::endl;
            if(DEBUG) std::cout << " numCurves: "   << numCurves   << std::endl;
            if(DEBUG) std::cout << " numSurfaces: " << numSurfaces << std::endl;
            if(DEBUG) std::cout << " numVolumes: "  << numVolumes  << std::endl;
            std::getline(buffer_Ent, line); // if ok, go to next line ..

            // Read Nodes
            if (numPoints>0) {
                for (size_t i=0; i<numPoints; i++)
                {
                    std::getline(buffer_Ent, line);
                    auto [ID, phys_ID] = get_entity(line,1); //1:node
                    std::cout << ID << phys_ID << std::endl;
                    //point2phys[ID] = phys_ID;
                }
            }
            // Read Curves
            if (numCurves>0) {
                for (size_t i=0; i<numCurves; i++)
                {
                    std::getline(buffer_Ent, line);
                    auto [ID, phys_ID] = get_entity(line,2); //2:curve
                    curve2phys[ID] = phys_ID;
                }
            }
            // Read Surfaces
            if (numSurfaces>0) {
                for (size_t i=0; i<numSurfaces; i++)
                {
                    std::getline(buffer_Ent, line);
                    auto [ID, phys_ID] = get_entity(line,3); //3:surface
                    surf2phys[ID] = phys_ID;
                }
            }
            // Read Volumes
            if (numVolumes>0) {
                for (size_t i=0; i<numVolumes; i++)
                {
                    std::getline(buffer_Ent, line);
                    auto [ID, phys_ID] = get_entity(line,4); //4:volume
                    volm2phys[ID] = phys_ID;
                }
            }
            // Clear Buffer
            buffer_Ent.str(std::string());

        } else {
            /**************************/
            // Read Partitioned Entities
            /**************************/
            std::stringstream buffer_PEnt;
            buffer_PEnt << PartEntities;
            
            size_t numPartitions = 0;
            buffer_PEnt >> numPartitions;

            std::getline(buffer_PEnt, line);
            size_t numGhostEntities = 0;
            buffer_PEnt >> numGhostEntities; // not used for the moment

            std::getline(buffer_PEnt, line);
            size_t numPoints  = 0;
            size_t numCurves  = 0;
            size_t numSurfaces= 0;
            size_t numVolumes = 0;
            buffer_PEnt >> numPoints >> numCurves >> numSurfaces >> numVolumes;

            if(DEBUG) std::cout << " numPartitions: " << numPartitions << std::endl;
            if(DEBUG) std::cout << " numPoints: "     << numPoints     << std::endl;
            if(DEBUG) std::cout << " numCurves: "     << numCurves     << std::endl;
            if(DEBUG) std::cout << " numSurfaces: "   << numSurfaces   << std::endl;
            if(DEBUG) std::cout << " numVolumes: "    << numVolumes    << std::endl;
            std::getline(buffer_PEnt, line); // if ok, go to next line ..

            // Read Nodes
            if (numPoints>0) {
                for (size_t i=0; i<numPoints; i++)
                {
                    std::getline(buffer_PEnt, line);
                    auto [chld_ID, prnt_ID, part_ID, phys_ID] = get_partitionedEntity(line,1); //1:node
                    point2part[chld_ID] = part_ID;
                    point2phys[chld_ID] = phys_ID;
                    point2geom[chld_ID] = prnt_ID;
                }
            }
            // Read Curves
            if (numCurves>0) {
                for (size_t i=0; i<numCurves; i++)
                {
                    std::getline(buffer_PEnt, line);
                    auto [chld_ID, prnt_ID, part_ID, phys_ID] = get_partitionedEntity(line,2); //2:curve
                    curve2part[chld_ID] = part_ID;
                    curve2phys[chld_ID] = phys_ID;
                    curve2geom[chld_ID] = prnt_ID;
                }
            }
            // Read Surfaces
            if (numSurfaces>0) {
                for (size_t i=0; i<numSurfaces; i++)
                {
                    std::getline(buffer_PEnt, line);
                    auto [chld_ID, prnt_ID, part_ID, phys_ID] = get_partitionedEntity(line,3); //3:surface
                    surf2part[chld_ID] = part_ID;
                    surf2phys[chld_ID] = phys_ID;
                    surf2geom[chld_ID] = prnt_ID;
                }
            }
            // Read Volumes
            if (numVolumes>0) {
                for (size_t i=0; i<numVolumes; i++)
                {
                    std::getline(buffer_PEnt, line);
                    auto [chld_ID, prnt_ID, part_ID, phys_ID] = get_partitionedEntity(line,4); //4:volume
                    volm2part[chld_ID] = part_ID;
                    volm2phys[chld_ID] = phys_ID;
                    volm2geom[chld_ID] = prnt_ID;
                }
            }
            // Clear Buffer
            buffer_PEnt.str(std::string());
        }

        /**************************/
        // Read Nodes
        /**************************/
        std::stringstream buffer_N;
        buffer_N << Nodes;
        
        size_t numNodesBlocks = 0;
        size_t numNodes       = 0;
        size_t minNodeIndex   = 0;
        size_t maxNodeIndex   = 0;
        buffer_N >> numNodesBlocks >> numNodes >> minNodeIndex >> maxNodeIndex;

        if(DEBUG) std::cout << " numNodesBlocks: " << numNodesBlocks   << std::endl;
        if(DEBUG) std::cout << " numNodes: "       << numNodes         << std::endl;
        if(DEBUG) std::cout << " minNodeIndex: "   << minNodeIndex-one << std::endl;
        if(DEBUG) std::cout << " maxNodeIndex: "   << maxNodeIndex-one << std::endl;

        //V = get_nodes(buffer_N,runNodesBlocks,numNodes);
        MArray<double,2> V = get_nodes(Nodes,numNodesBlocks,numNodes);

        // Report the number for nodes found
        std::cout << "Total vertices found = " << numNodes << std::endl;

        // Clear Buffer
        buffer_N.str(std::string());

        /**************************/
        // Read Elements
        /**************************/
        std::stringstream buffer_E;
        buffer_E << Elements;

        size_t numEntBlocks = 0;
        size_t numElements  = 0;
        size_t minElemIndex = 0; // not needed
        size_t maxElemIndex = 0; // not needed
        buffer_E >> numEntBlocks >> numElements >> minElemIndex >> maxElemIndex;

        if(DEBUG) std::cout << " numEntBlocks: " << numEntBlocks     << std::endl;
        if(DEBUG) std::cout << " numElements: "  << numElements      << std::endl;
        if(DEBUG) std::cout << " minElemIndex: " << minElemIndex-one << std::endl;
        if(DEBUG) std::cout << " maxElemIndex: " << maxElemIndex-one << std::endl;
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
        size_t e15 = 0; // point Element counter
        
        // Read elements block
        for (size_t Ent=0; Ent<numEntBlocks; Ent++)
        {
            // Read Block parameters
            std::getline(buffer_E, line);

            std::stringstream s_stream(line);
            size_t entityDim; // 0:point, 1:curve, 2:surface, 3:volume
            size_t entityTag; // this is: Entity.ID | Entity.child_ID
            size_t elementType; // 1:line, 2:triangle, 4:tetrahedron, 15:point
            size_t numElementsInBlock;
            s_stream >> entityDim >> entityTag >> elementType >> numElementsInBlock;

            // Read Elements in block
            for (size_t i=0; i<numElementsInBlock; i++)
            {
                std::getline(buffer_E, line);
                auto line_data = str2size_t(line);
                switch (elementType) // <-- Should use entityDim, but we only search 4 type of elements
                {
                case 1: /* Line elements */
                    e1 = e1 + 1; // update element counter
                    break;
                case 2: /* triangle elements */
                    e2 = e2 + 1; // update element counter
                    break;
                case 4: /* tetrahedron elements */
                    e4 = e4 + 1; // update element counter
                    break;
                case 15: /* point elements */
                    e15 = e15 + 1; // update element counter
                    break;
                default:
                    break;
                }
            }
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
        
    } else {
        std::cout << "Could not open file: " << mesh_file << std::endl; 
        std::exit(-1);
    }
    // If everything goes well ...
    return 0;
}

#endif