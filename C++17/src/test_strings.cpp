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

    // Print comments 
    bool DEBUG = false; // print all all warnings and stage messages

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

        // DEBUG
        if(DEBUG) std::cout << MeshFormat << std::endl;
        if(DEBUG) std::cout << PhysicalNames << std::endl;
        if(DEBUG) std::cout << Entities << std::endl;
        if(DEBUG) std::cout << PartEntities << std::endl;
        if(DEBUG) std::cout << Nodes << std::endl;
        if(DEBUG) std::cout << Elements << std::endl;

        // Read the first line from Mesh format
        buffer.str(Nodes);
        std::string line;
        for(int i=0; i<4; i++){
            std::getline(buffer,line);
            std::cout << line << std::endl;
        }

        // Clear buffer-variable
        buffer.str(std::string());  // this is equivalent to: buffer.str("");
    } else {
        //std::cout << "Could not open file: " << mesh_file << std::endl; 
        std::exit(-1);
    }
    // If everything goes well ...
    return 0;
}