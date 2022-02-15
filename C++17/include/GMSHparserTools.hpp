#ifndef GMSH_PARSER_TOOLS
#define GMSH_PARSER_TOOLS

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//   Common tools required for reading GMSH-file in format v2.2 and v4.1 
//
//      Coded by Manuel A. Diaz @ Pprime | Univ-Poitiers, 2022.01.21
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Globals.hpp"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Get a sub-string from the main string-buffer using two unique delimiters:
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Get a sub-vector from a std::vector using two delimiters:
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<size_t> extractVectorBetween(
    const std::vector<size_t> Vec,
    const size_t first_index,
    const size_t last_index)
{
    std::vector<size_t>::const_iterator first = Vec.begin() + first_index;
    std::vector<size_t>::const_iterator last = Vec.begin() + last_index;
    std::vector<size_t> subVec(first, last);

    return subVec;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Build convention map of Boundary conditions for ParadigmS:
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Element structure to be acquired:
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct element
{
    std::vector<size_t> EToV;
    std::vector<int> phys_tag;
    std::vector<int> geom_tag;
    std::vector<int> part_tag;
    std::vector<int> Etype;
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Split and input strings and return a vector with all double values
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<double> str2double(const std::string &str)
{
    std::stringstream ss(str);
    std::istream_iterator<double> begin(ss);
    std::istream_iterator<double> end;
    std::vector<double> values(begin,end);
    return values;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Split and input strings and return a vector with all size_t values
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<size_t> str2size_t(const std::string &str)
{
    std::stringstream ss(str);
    std::istream_iterator<size_t> begin(ss);
    std::istream_iterator<size_t> end;
    std::vector<size_t> values(begin,end);
    return values;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Get entity Tag and its associated physical Tag.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Get partitioned entity Tags and its associated physical Tag.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                physicalTag = int(vector[11+numPartitionsTags]);
            }
            break;
    }
    return {entityTag, parentTag, partitionTag, physicalTag};
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Get nodes from block system (GMSG format 4.1)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MArray<double,2> get_nodes(const std::string &Nodes, const size_t &numNodesBlocks, const size_t &numNodes, const size_t Dim) 
{
    // Allocate output
    MArray<double,2> V({numNodes,Dim},0); // [x(:),y(:),z(:)]

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
            if (Dim==2) {stream >> V(n,0) >> V(n,1);}
            if (Dim==3) {stream >> V(n,0) >> V(n,1) >> V(n,2);}
            n = n+1; // Update node counter
        }
        // Delete temporary new-arrays
        delete [] nodeTag;
    }
    return V;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Get nodes from block system (GMSH format 2.2)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MArray<double,2> get_nodes(const std::string &Nodes, const size_t &numNodes, const size_t Dim) 
{
    // Allocate output
    MArray<double,2> V({numNodes,Dim},0); // [x(:),y(:),z(:)]

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
        if (Dim==2) {stream >> nodeTag >> V(i,0) >> V(i,1);}
        if (Dim==3) {stream >> nodeTag >> V(i,0) >> V(i,1) >> V(i,2);}
    }
    return V;
}
#endif