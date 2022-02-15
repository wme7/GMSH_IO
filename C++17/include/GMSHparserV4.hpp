#ifndef GMSH_PARSER_V4
#define GMSH_PARSER_V4

// local libs
#include "Globals.hpp"
#include "GMSHparserTools.hpp"

int GMSHparserV4(std::string mesh_file) 
{
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    //     Extract entities contained in a single GMSH-file in format v4.1 
    //
    //      Coded by Manuel A. Diaz @ Pprime | Univ-Poitiers, 2022.01.21
    //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    // Example call: GMSHparserV4('filename.msh')
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
    // -- READ GMSH V4.1 FORMAT --
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

        if (single_domain) {
            /**************************/
            // 3. Read Entities
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
                    point2phys[ID] = phys_ID;
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
            /*******************************/
            // 4. Read Partitioned Entities
            /*******************************/
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
        // 5. Read Nodes
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
        // 6. Read Elements
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
        size_t numE1 = 0; // Lines Element counter
        size_t numE2 = 0; // Triangle Element counter
        size_t numE4 = 0; // Tetrahedron Element counter
        size_t numE15= 0; // Point Element counter
        
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
                std::getline(buffer_E,line);
                auto line_data = str2size_t(line);
                switch (elementType) // <-- Should use entityDim, but we only search 4 type of elements
                {
                case 1: /* Line elements */
                    numE1 = numE1 + 1; // update element counter
                    LE.Etype.push_back(elementType);
                    LE.EToV.push_back(line_data[1]-one);
                    LE.EToV.push_back(line_data[2]-one);
                    LE.phys_tag.push_back(curve2phys[entityTag]);
                    if (not(single_domain)) {
                        LE.geom_tag.push_back(curve2geom[entityTag]);
                        LE.part_tag.push_back(curve2part[entityTag]);
                    } else {
                        LE.geom_tag.push_back(entityTag);
                    }
                    break;
                case 2: /* triangle elements */
                    numE2 = numE2 + 1; // update element counter
                    SE.Etype.push_back(elementType);
                    SE.EToV.push_back(line_data[1]-one);
                    SE.EToV.push_back(line_data[2]-one);
                    SE.EToV.push_back(line_data[3]-one);
                    SE.phys_tag.push_back(surf2phys[entityTag]);
                    if (not(single_domain)) {
                        SE.geom_tag.push_back(surf2geom[entityTag]);
                        SE.part_tag.push_back(surf2part[entityTag]);
                    } else {
                        SE.geom_tag.push_back(entityTag);
                    }
                    break;
                case 4: /* tetrahedron elements */
                    numE4 = numE4 + 1; // update element counter
                    VE.Etype.push_back(elementType);
                    VE.EToV.push_back(line_data[1]-one);
                    VE.EToV.push_back(line_data[2]-one);
                    VE.EToV.push_back(line_data[3]-one);
                    VE.EToV.push_back(line_data[4]-one);
                    VE.phys_tag.push_back(volm2phys[entityTag]);
                    if (not(single_domain)) {
                        VE.geom_tag.push_back(volm2geom[entityTag]);
                        VE.part_tag.push_back(volm2part[entityTag]);
                    } else {
                        VE.geom_tag.push_back(entityTag);
                    }
                    break;
                case 15: /* point elements */
                    numE15 = numE15 + 1; // update element counter
                    PE.Etype.push_back(elementType);
                    PE.EToV.push_back(line_data[1]-one);
                    PE.phys_tag.push_back(point2phys[entityTag]);
                    if (not(single_domain)) {
                        PE.geom_tag.push_back(point2geom[entityTag]);
                        PE.part_tag.push_back(point2part[entityTag]);
                    } else {
                        PE.geom_tag.push_back(entityTag);
                    }
                    break;
                default:
                    std::cout << "ERROR: element type not in list"<< std::endl;
                    std::exit(-1);
                    break;
                }
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
        // 7. Parse FEM data
        /**************************/

        // Save to MdimArrays
        MArray<size_t,2> PEToV({numE15,1},PE.EToV); //PEToV.print();
        MArray<size_t,2> LEToV({ numE1,2},LE.EToV); //LEToV.print();
        MArray<size_t,2> SEToV({ numE2,3},SE.EToV); //SEToV.print();
        MArray<size_t,2> VEToV({ numE4,4},VE.EToV); //VEToV.print();

        // Save to Numpy Array
        cnpy::npz_save("../../Python3/FEMmeshV4.npy",  "V"  , V.data(), {numNodes,3},"w");
        cnpy::npz_save("../../Python3/FEMmeshV4.npy","PEToV",PEToV.data(),{numE15,1},"a");
        cnpy::npz_save("../../Python3/FEMmeshV4.npy","LEToV",LEToV.data(),{numE1, 2},"a");
        cnpy::npz_save("../../Python3/FEMmeshV4.npy","SEToV",SEToV.data(),{numE2, 3},"a");
        cnpy::npz_save("../../Python3/FEMmeshV4.npy","VEToV",VEToV.data(),{numE4, 4},"a");
        
    } else {
        std::cout << "ERROR: Could not open file: " << mesh_file << std::endl; 
        std::exit(-1);
    }
    // If everything goes well ...
    return 0;
}
#endif