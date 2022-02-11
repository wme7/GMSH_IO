#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Get a sub-string from the main string-buffer using two unique delimiters:
std::string extractBetween(
    const std::string &buffer,
    const std::string &start_delimiter,
    const std::string &stop_delimiter)
{
    unsigned first_delim_pos = buffer.find(start_delimiter) + 1 ; // add '\n' char
    unsigned end_pos_of_first_delim = first_delim_pos + start_delimiter.length();
    unsigned last_delim_pos = buffer.find(stop_delimiter) - 1; // minus '\n' char
 
    return buffer.substr(end_pos_of_first_delim,
        last_delim_pos - end_pos_of_first_delim);
}
 
int main() {

    // We assume that the maximun dimension within the physical groups
    // will be the maximun dimension of the problem
    size_t phys_MAX_DIM = 0; // initial guess
    size_t phys_MIN_DIM = 3; // initial guess

    // Print comments 
    bool DEBUG = true; // print all all warnings and stage messages

    // Index correction factor
    size_t one = 1; // if we set one=0, one recovers the original indexes of GMSH.

    // Last by not least, open meshfile:    
    std::ifstream file("../../meshes/square_v4_2P.msh");
    //std::ifstream fichier(mesh_file); 

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
        std::string MeshFormat    = extractBetween(buffer.str(),"$MeshFormat","$EndMeshFormat");
        std::string PhysicalNames = extractBetween(buffer.str(),"$PhysicalNames","$EndPhysicalNames");
        std::string Entities      = extractBetween(buffer.str(),"$Entities","$EndEntities");
        std::string PartEntities  = extractBetween(buffer.str(),"$PartitionedEntities","$EndPartitionedEntities");
        std::string Nodes         = extractBetween(buffer.str(),"$Nodes","$EndNodes");
        std::string Elements      = extractBetween(buffer.str(),"$Elements","$EndElements");

        // Sanity check
        if (  MeshFormat.empty() ) {std::cout << " Error - Wrong File Format!" << std::endl; std::exit(-1);}
        if (PhysicalNames.empty()) {std::cout << " Error - No Physical names!" << std::endl; std::exit(-1);}
        if (   Entities.empty()  ) {std::cout << " Error - No Entities found!" << std::endl; std::exit(-1);}
        if (    Nodes.empty()    ) {std::cout << " Error - Nodes are missing!" << std::endl; std::exit(-1);}
        if (   Elements.empty()  ) {std::cout << " Error - No Elements found!" << std::endl; std::exit(-1);}

        //if(DEBUG) std::cout << MeshFormat << std::endl;
        //if(DEBUG) std::cout << PhysicalNames << std::endl;
        //if(DEBUG) std::cout << Entities << std::endl;
        //if(DEBUG) std::cout << PartEntities << std::endl;
        //if(DEBUG) std::cout << Nodes << std::endl;
        //if(DEBUG) std::cout << Elements << std::endl;

        // For reading the buffer line by line
        std::string line;

        /**************************/
        // 1. Read MeshFormat
        /**************************/
        std::stringstream buffer_MF;    
        buffer_MF << MeshFormat;

            double version = 1.0; size_t format=0, size=0;
            buffer_MF >> version >> format >> size;

            if(DEBUG) std::cout << "version: " << version << std::endl;
            if(DEBUG) std::cout << "format: " << format << std::endl;
            if(DEBUG) std::cout << "size: " << size << std::endl;

            if (version != 4.1)
            {
                std::cout << " Error - incorrect GMSH version" << std::endl; 
                std::exit(-1);
            }
            if (format != 0)
            {
                std::cout << " Error - not ASCII file format" << std::endl; 
                std::exit(-1);
            }

        /**************************/
        // Read PhysicalNames
        /**************************/
        std::stringstream buffer_PN;
        buffer_PN << PhysicalNames;

            // NOTE:
            // The lines in the PhysicalNames section should look like the following:
            // 
            // $PhysicalNames
            // 5
            // 2 1 "left BC"
            // 2 3 "right BC"
            // 2 4 "top BC"
            // 2 5 "bottom BC"
            // 3 2 "volume"
            // $EndPhysicalNames

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

                if(DEBUG) std::cout << " gmsh physical: " << phys_id << std::endl;
                if(DEBUG) std::cout << " dimension: " << phys_dim << std::endl;
                if(DEBUG) std::cout << " name: " << phys_name << std::endl;

                // Not sure if this is true for all Gmsh files, but
                // my test file has quotes around the phys_name string.
                // So let's erase any quotes now...
                phys_name.erase(std::remove(phys_name.begin(), phys_name.end(), '"'), phys_name.end());

                // // Record this ID for later assignment of subdomain/sideset names.
                // gmsh_physicals[phys_id] = std::make_pair(phys_dim, phys_name);

                // if(DEBUG) std::cout << " phys_dim: " << std::get<0>(gmsh_physicals[phys_id]) << std::endl;
                // if(DEBUG) std::cout << " phys_name: " << std::get<1>(gmsh_physicals[phys_id]) << std::endl;

                // Search for the maximun and minimun dimensions of the computational domain:
                phys_MAX_DIM = phys_MAX_DIM > phys_dim ? phys_MAX_DIM : phys_dim;
                phys_MIN_DIM = phys_MIN_DIM < phys_dim ? phys_MIN_DIM : phys_dim;

                if(DEBUG) std::cout << " MAX_DIM: " << phys_MAX_DIM << std::endl;
                if(DEBUG) std::cout << " MIN_DIM: " << phys_MIN_DIM << std::endl;
            }

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


        /**************************/
        // Read Partitioned Entities
        /**************************/
        if (not(PartEntities.empty()))
        {
            std::cout << "Partitioned domain detected" << std::endl;
            std::stringstream buffer_PEnt;
            buffer_PEnt << PartEntities;
            
                size_t numPartitions = 0;
                buffer_PEnt >> numPartitions;

                std::getline(buffer_PEnt, line);
                size_t numGhostEntities = 0;
                buffer_PEnt >> numGhostEntities; // not used for the moment

                std::getline(buffer_PEnt, line);
                size_t numPartPoints  = 0;
                size_t numPartCurves  = 0;
                size_t numPartSurfaces= 0;
                size_t numPartVolumes = 0;
                buffer_PEnt >> numPartPoints >> numPartCurves >> numPartSurfaces >> numPartVolumes;

                if(DEBUG) std::cout << " numPartitions: " << numPartitions   << std::endl;
                if(DEBUG) std::cout << " numPoints: "     << numPartPoints   << std::endl;
                if(DEBUG) std::cout << " numCurves: "     << numPartCurves   << std::endl;
                if(DEBUG) std::cout << " numSurfaces: "   << numPartSurfaces << std::endl;
                if(DEBUG) std::cout << " numVolumes: "    << numPartVolumes  << std::endl;

                // Clear Buffer
                buffer_PEnt.str(std::string());
        } else {
            // Do nothing !
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

        /**************************/
        // Read Elements
        /**************************/
        std::stringstream buffer_E;
        buffer_E << Elements;

            size_t numElemBlocks= 0;
            size_t numElements  = 0;
            size_t minElemIndex = 0;
            size_t maxElemIndex = 0;
            buffer_E >> numElemBlocks >> numElements >> minElemIndex >> maxElemIndex;

            if(DEBUG) std::cout << " numElemBlocks: " << numElemBlocks    << std::endl;
            if(DEBUG) std::cout << " numElements: "   << numElements      << std::endl;
            if(DEBUG) std::cout << " minElemIndex: "  << minElemIndex-one << std::endl;
            if(DEBUG) std::cout << " maxElemIndex: "  << maxElemIndex-one << std::endl;
        
        /**************************/
        // Clear Buffers
        /**************************/
        buffer_MF.str(std::string());  // this is equivalent to: buffer.str("");
        buffer_PN.str(std::string());
        buffer_Ent.str(std::string());
        buffer_N.str(std::string());
        buffer_E.str(std::string());
    } else {
        //std::cout << "Could not open file: " << mesh_file << std::endl; 
        std::exit(-1);
    }
    // If everything goes well ...
    return 0;
}