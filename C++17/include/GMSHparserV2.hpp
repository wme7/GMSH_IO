#ifndef GMSH_PARSER_V2
#define GMSH_PARSER_V2

// local libs
#include "Globals.hpp"
#include "GMSHparserTools.hpp"

int GMSHparserV2(std::string mesh_file) 
{
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    //     Extract entities contained in a single GMSH-file in format v2.2 
    //
    //      Coded by Manuel A. Diaz @ Pprime | Univ-Poitiers, 2022.01.21
    //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    // Example call: GMSHparserV2('filename.msh')
    //
    // Output:
    //     V: the vertices (nodes coordinates) -- (Nx3) array
    //    VE: volumetric elements (tetrahedrons) -- structure
    //    SE: surface elements (triangles,quads) -- structure
    //    LE: curvilinear elements (lines/edges) -- structure
    //    PE: point elements (singular vertices) -- structure
    //    mapPhysNames: maps phys.tag --> phys.name  -- map structure
    //    info: version, format, endian test -- structure
    //
    // Note: This parser is designed to capture the following elements:
    //
    // Point:                Line:                   Triangle:
    //
    //        v                                              v
    //        ^                                              ^
    //        |                       v                      |
    //        +----> u                ^                      2
    //       0                        |                      |`\.
    //                                |                      |  `\.
    //                          0-----+-----1 --> u          |    `\.
    //                                                       |      `\.
    // Tetrahedron:                                          |        `\.
    //                                                       0----------1 --> u
    //                    v
    //                   ,
    //                  /
    //               2
    //             ,/|`\                    Based on the GMSH guide 4.9.4
    //           ,/  |  `\                  This are lower-order elements 
    //         ,/    '.   `\                identified as:
    //       ,/       |     `\                  E-1 : 2-node Line 
    //     ,/         |       `\                E-2 : 3-node Triangle
    //    0-----------'.--------1 --> u         E-4 : 4-node tetrahedron
    //     `\.         |      ,/                E-15: 1-node point
    //        `\.      |    ,/
    //           `\.   '. ,/                Other elements can be added to 
    //              `\. |/                  this parser by modifying the 
    //                 `3                   Read Elements stage.
    //                    `\.
    //                       ` w            Happy coding ! M.D. 02/2022.
    //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //-------------------------------------
    // - Global Parameters of the parser -
    //-------------------------------------

    // Nominal dimension of the meshed geometry
    size_t phys_DIM = 0; // initial guess

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
        std::string Nodes         = extractBetween(buffer.str(),"$Nodes\n","\n$EndNodes");
        std::string Elements      = extractBetween(buffer.str(),"$Elements\n","\n$EndElements");

        // Sanity check
        if (  MeshFormat.empty() ) {std::cout << " Error - Wrong File Format!" << std::endl; std::exit(-1);}
        if (PhysicalNames.empty()) {std::cout << " Error - No Physical names!" << std::endl; std::exit(-1);}
        if (    Nodes.empty()    ) {std::cout << " Error - Nodes are missing!" << std::endl; std::exit(-1);}
        if (   Elements.empty()  ) {std::cout << " Error - No Elements found!" << std::endl; std::exit(-1);}

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
        MArray<double,2> V = get_nodes(Nodes,numNodes,phys_DIM); //V.print();

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
        size_t numE1 = 0; // Lines Element counter
        size_t numE2 = 0; // Triangle Element counter
        size_t numE4 = 0; // Tetrahedron Element counter
        size_t numE15= 0; // Point Element counter

        // Read elements serially
        for (size_t i=0; i<numElements; i++) 
        {
            // Read Block parameters
            std::getline(buffer_E, line);
            auto n = str2size_t(line);
            //size_t elementID = n(0); % we use a local numbering instead
            size_t elementType = n[1];
            size_t numberOfTags = n[2];
            switch (elementType)
            {
            case 1: // Line elements
                numE1 = numE1 + 1; // update element counter
                LE.Etype.push_back(elementType);
                LE.EToV.push_back(n[2+numberOfTags+1]-one);
                LE.EToV.push_back(n[2+numberOfTags+2]-one);
                if (numberOfTags>0) {
                    auto tags = extractVectorBetween(n,3,3+numberOfTags);
                    // tags(1) : physical entity to which the element belongs
                    // tags(2) : elementary number to which the element belongs
                    // tags(3) : number of partitions to which the element belongs
                    // tags(4) : partition id number
                    if (tags.size()>=1) {
                        LE.phys_tag.push_back(tags[0]);
                        if (tags.size()>=2) {
                            LE.geom_tag.push_back(tags[1]);
                            if (tags.size()>=4) {
                                LE.part_tag.push_back(tags[3]);
                            }
                        }
                    }
                }
                break;
            case 2: // Triangle elements
                numE2 = numE2 + 1; // update element counter
                SE.Etype.push_back(elementType);
                SE.EToV.push_back(n[2+numberOfTags+1]-one);
                SE.EToV.push_back(n[2+numberOfTags+2]-one);
                SE.EToV.push_back(n[2+numberOfTags+3]-one);
                if (numberOfTags>0) {
                    auto tags = extractVectorBetween(n,3,3+numberOfTags);
                    // tags(1) : physical entity to which the element belongs
                    // tags(2) : elementary number to which the element belongs
                    // tags(3) : number of partitions to which the element belongs
                    // tags(4) : partition id number
                    if (tags.size()>=1) {
                        SE.phys_tag.push_back(tags[0]);
                        if (tags.size()>=2) {
                            SE.geom_tag.push_back(tags[1]);
                            if (tags.size()>=4) {
                                SE.part_tag.push_back(tags[3]);
                            }
                        }
                    }
                }
                break;
            case 4: // Tetrahedron elements
                numE4 = numE4 + 1; // update element counter
                VE.Etype.push_back(elementType);
                VE.EToV.push_back(n[2+numberOfTags+1]-one);
                VE.EToV.push_back(n[2+numberOfTags+2]-one);
                VE.EToV.push_back(n[2+numberOfTags+3]-one);
                VE.EToV.push_back(n[2+numberOfTags+4]-one);
                if (numberOfTags>0) {
                    auto tags = extractVectorBetween(n,3,3+numberOfTags);
                    // tags(1) : physical entity to which the element belongs
                    // tags(2) : elementary number to which the element belongs
                    // tags(3) : number of partitions to which the element belongs
                    // tags(4) : partition id number
                    if (tags.size()>=1) {
                        VE.phys_tag.push_back(tags[0]);
                        if (tags.size()>=2) {
                            VE.geom_tag.push_back(tags[1]);
                            if (tags.size()>=4) {
                                VE.part_tag.push_back(tags[3]);
                            }
                        }
                    }
                }
                break;
            case 15: // Point elements
                numE15 = numE15 + 1; // update element counter
                PE.Etype.push_back(elementType);
                PE.EToV.push_back(n[2+numberOfTags+1]-one);
                if (numberOfTags>0) {
                    auto tags = extractVectorBetween(n,3,3+numberOfTags);
                    // tags(1) : physical entity to which the element belongs
                    // tags(2) : elementary number to which the element belongs
                    // tags(3) : number of partitions to which the element belongs
                    // tags(4) : partition id number
                    if (tags.size()>=1) {
                        PE.phys_tag.push_back(tags[0]);
                        if (tags.size()>=2) {
                            PE.geom_tag.push_back(tags[1]);
                            if (tags.size()>=4) {
                                PE.part_tag.push_back(tags[3]);
                            }
                        }
                    }
                }
                break;
            default:
                std::cout << "ERROR: element type not in list"<< std::endl;
                std::exit(-1);
                break;
            }
        }
        std::cout << "Total point-elements found = " << numE15 << std::endl;
        std::cout << "Total curve-elements found = " << numE1 << std::endl;
        std::cout << "Total surface-elements found = " << numE2 << std::endl;   
        std::cout << "Total volume-elements found = " << numE4 << std::endl;
        // Sanity check
        if (numElements != (numE15+numE1+numE2+numE4)) {
            std::cout << "Total number of elements missmatch!"<< std::endl;
            std::exit(-1);
        }
        // Clear Buffer
        buffer_E.str(std::string());

        /**************************/
        // 5. Parse FEM data
        /**************************/

        // Save to MdimArrays
        MArray<size_t,2> PEToV({numE15,1},PE.EToV); //PEToV.print();
        MArray<size_t,2> LEToV({ numE1,2},LE.EToV); //LEToV.print();
        MArray<size_t,2> SEToV({ numE2,3},SE.EToV); //SEToV.print();
        MArray<size_t,2> VEToV({ numE4,4},VE.EToV); //VEToV.print();

        std::cout << "LE.phys_tag\n" << std::endl;
        std::copy(LE.phys_tag.begin(), LE.phys_tag.end(), std::ostream_iterator<int>(std::cout,"\n"));
        std::cout << "LE.geom_tag\n" << std::endl;
        std::copy(LE.geom_tag.begin(), LE.geom_tag.end(), std::ostream_iterator<int>(std::cout,"\n"));
        std::cout << "LE.part_tag\n" << std::endl;
        std::copy(LE.part_tag.begin(), LE.part_tag.end(), std::ostream_iterator<int>(std::cout,"\n"));

        // Save to Numpy Array
        cnpy::npz_save("../../Python3/FEMmeshV2.npz","V",V.data(),{numNodes,phys_DIM},"w");
        cnpy::npz_save("../../Python3/FEMmeshV2.npz","PEToV",PEToV.data(),{numE15,1},"a");
        cnpy::npz_save("../../Python3/FEMmeshV2.npz","LEToV",LEToV.data(),{numE1, 2},"a");
        cnpy::npz_save("../../Python3/FEMmeshV2.npz","SEToV",SEToV.data(),{numE2, 3},"a");
        cnpy::npz_save("../../Python3/FEMmeshV2.npz","VEToV",VEToV.data(),{numE4, 4},"a");

    } else {
        std::cout << "ERROR: Could not open file: " << mesh_file << std::endl; 
        std::exit(-1);
    }
    // If everything goes well ...
    return 0;
}
#endif