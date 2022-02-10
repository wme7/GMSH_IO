#ifndef GMSH_PARSER
#define GMSH_PARSER

// local libs
#include "Globals.hpp"
#include "Node.hpp"
#include "Point.hpp"
#include "Line.hpp"
#include "Triangle.hpp"
#include "Tetrahedron.hpp"

// The GMSHparser functions is here used to read 2D/3D unstructured meshes in a uniform manner.
// It provides access to properties of the mesh such as number of elements and their vertex. 
// It also provides the following standarized arrays: V, EToVc, EToV, BEToV, EToP, BEToP.
//
// Coded by MD in UPMC 26.04.2019

FEMmesh GMSHparser(std::string mesh_file)
{
    //-------------------------------------
    //  0 - INITIALIZATION
    //-------------------------------------
    // This is a serial-only process for now;
    // the Mesh should be read on processor 0 and
    // broadcast later
    
    FEMmesh m; // See more details in Globals

    // Mapping from physical id -> (physical dim, physical name) pairs.
    // These can refer to either "sidesets" or "subdomains"; we need to
    // wait until the Mesh has been read to know which is which.  
    std::map<size_t, std::pair<size_t, std::string>> gmsh_physicals;
    std::map<size_t, std::string> gmsh_physical_names;

    // We assume that the maximun dimension within the physical groups
    // will be the maximun dimension of the problem
    size_t phys_MAX_DIM = 0; // initial guess
    size_t phys_MIN_DIM = 3; // initial guess

    // Initilize Nodes and Element containers
    std::vector<Node<double>> Nodes;
    std::vector<Point> Points;
    std::vector<Line> Lines;
    std::vector<Triangle> Faces;
    std::vector<Line> BoundaryElements2d;
    std::vector<Triangle> Elements2d;
    std::vector<Triangle> BoundaryElements3d;
    std::vector<Tetrahedron> Elements3d;

    // Print comments 
    bool DEBUG = false; // print all all warnings and stage messages

    // Index correction factor
    size_t one = 1; // if we set one=0, one recovers the original indexes of GMSH.

    // Last by not least, open meshfile:
    std::ifstream fichier(mesh_file); 

    //-------------------------------------
    //  1 - READ GMSH FORMAT
    //-------------------------------------
    if (fichier.is_open())
    {
        if(DEBUG) std::cout << "Starting to read mesh in file: "<< mesh_file << std::endl;
        
        // For reading the file line by line
        std::string ligne;

        while (true)
        {
            std::getline(fichier, ligne); // get whatever is in the line
            // Process line
            if (ligne.find("$MeshFormat") == static_cast<std::string::size_type>(0))
            {
                double version = 1.0; size_t format=0, size=0;
                fichier >> version >> format >> size;

                if(DEBUG) std::cout << "version: " << version << std::endl;
                if(DEBUG) std::cout << "format: " << format << std::endl;
                if(DEBUG) std::cout << "size: " << size << std::endl;
                
                if (version != 2.2)
                {
                    std::cout << " Error - not supported msh file version" << std::endl; std::exit(-1);
                }
                if (format != 0)
                {
                    std::cout << " Error - not ASCII file format" << std::endl; std::exit(-1);
                }
            }
            // Read the "PhysicalNames" section
            else if (ligne.find("$PhysicalNames") == static_cast<std::string::size_type>(0))
            {
                // Read in the number of physical groups to expect in the file.
                size_t num_physical_groups = 0;
                fichier >> num_physical_groups;

                if(DEBUG) std::cout << "num_physical_groups: " << num_physical_groups << std::endl;

                // Read rest of line including newline character.
                std::getline(fichier, ligne);

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

                for (size_t i=0; i<num_physical_groups; ++i)
                {
                    // Read an entire line of the PhysicalNames section.
                    std::getline(fichier, ligne);

                    // Use an istringstream to extract the physical
                    // dimension, physical id, and physical name from
                    // this line.
                    std::istringstream s_stream(ligne);
                    size_t phys_dim, phys_id;
                    std::string phys_name;

                    s_stream >> phys_dim >> phys_id >> phys_name;

                    if(DEBUG) std::cout << " gmsh physical: " << phys_id << std::endl;
                    if(DEBUG) std::cout << " dimension: " << phys_dim << std::endl;
                    if(DEBUG) std::cout << " name: " << phys_name << std::endl;

                    // Not sure if this is true for all Gmsh files, but
                    // my test file has quotes around the phys_name string.
                    // So let's erase any quotes now...
                    phys_name.erase(std::remove(phys_name.begin(), phys_name.end(), '"'), phys_name.end());

                    // Record this ID for later assignment of subdomain/sideset names.
                    gmsh_physicals[phys_id] = std::make_pair(phys_dim, phys_name);

                    if(DEBUG) std::cout << " phys_dim: " << std::get<0>(gmsh_physicals[phys_id]) << std::endl;
                    if(DEBUG) std::cout << " phys_name: " << std::get<1>(gmsh_physicals[phys_id]) << std::endl;

                    // Search for the maximun and minimun dimensions of the computational domain:
                    phys_MAX_DIM = phys_MAX_DIM > phys_dim ? phys_MAX_DIM : phys_dim;
                    phys_MIN_DIM = phys_MIN_DIM < phys_dim ? phys_MIN_DIM : phys_dim;

                    if(DEBUG) std::cout << "MAX_DIM: " << phys_MAX_DIM << std::endl;
                    if(DEBUG) std::cout << "MIN_DIM: " << phys_MIN_DIM << std::endl;
                }
                // Verify that is the actual end of the section
                std::getline(fichier, ligne);
                if(ligne.find("$EndPhysicalNames") != static_cast<std::string::size_type>(0))
                {
                    std::cout << " Error - Something is wrong with the mesh" << std::endl; std::exit(-1);
                } 
            }
            // Read the node block
            else if (ligne.find("$Nodes") == static_cast<std::string::size_type>(0))
            {
                // Read in the number of nodes to expect in the file.
                fichier >> m.num_nodes;

                if(DEBUG) std::cout << "num_nodes: " << m.num_nodes << std::endl;

                // read in the nodal coordinates and form points.
                size_t id; 
                double x, y, z;

                // Read Nodes and add them to the container
                for (size_t i=0; i<m.num_nodes; ++i)
                {
                    fichier >> id >> x >> y >> z;
                    Node<double> node(id-1,x,y,z); 
                    Nodes.push_back(node);
                    if(DEBUG) std::cout << id-one << " " << x << " " << y << " " << z << std::endl;
                }
                // get the end of this line
                std::getline(fichier, ligne);

                // Verify that this is the actual end of the section
                std::getline(fichier, ligne);
                if(ligne.find("$EndNodes") != static_cast<std::string::size_type>(0))
                {
                    std::cout << " Error - Something is wrong with the nodes" << std::endl; std::exit(-1);
                } 
            }
            // Read the element block
            else if (ligne.find("$Elements") == static_cast<std::string::size_type>(0))
            {
                // read how many elements are there, and reserve space in the mesh
                size_t num_elements = 0;
                fichier >> num_elements;

                if(DEBUG) std::cout << "num_elements: " << num_elements << std::endl;

                // finish reading the current line
                std::getline(fichier, ligne);

                // Elements2d properties.
                size_t id, elementType; int nTags, negativeTag;
                size_t physTag, geometryTag, num_PartitionTag=1, partitionTag=1;

                // store the partition Tag for each element/control object
                std::vector<size_t> K_partitionTags; 
                std::vector<size_t> KB_partitionTags;
                std::vector<size_t> CO_Points_partitionTags;
                std::vector<size_t> CO_Lines_partitionTags; 
                std::vector<size_t> CO_Faces_partitionTags;  

                // Read elements depending on the type of mesh detected:
                if (phys_MAX_DIM==2) // We assume this is a 2-dimensional simulation problem.
                {
                    m.spatialDimension = 2;
                    m.Nfaces = 3;
                    if(DEBUG) std::cout << "Domain for 2-d test has been found." << std::endl;

                    // Read all types of elements and add them to their respective container
                    for (size_t i=0; i<num_elements; ++i)
                    {
                        fichier >> id >> elementType >> nTags;
                        if(DEBUG) std::cout << id-one << " " << elementType << " " << nTags << " ";

                        if(nTags==2)  // single domain
                        {  
                            fichier >> physTag >> geometryTag; num_PartitionTag=1; partitionTag=1;
                            if(DEBUG) std::cout << physTag << " " << geometryTag << " ";
                        } 
                        if(nTags==4) // partitioned domain
                        {   
                            fichier >> physTag >> geometryTag >> num_PartitionTag >> partitionTag;
                            if(DEBUG) std::cout << physTag << " " << geometryTag << " " << num_PartitionTag << " " << partitionTag-one << " ";
                        }
                        if(nTags>4) // partitioned domain
                        {   
                            fichier >> physTag >> geometryTag >> num_PartitionTag >> partitionTag;
                            if(DEBUG) std::cout << physTag << " " << geometryTag << " " << num_PartitionTag << " " << partitionTag-one << " ";

                            for (int i=0; i<(nTags-4); i++)
                            {
                                fichier >> negativeTag;
                                if(DEBUG) std::cout << negativeTag << " ";
                            }
                        }
                        if(elementType==1) // Edge (Boundary) element
                        {   
                            // Nodes and boundaries type 
                            size_t node1, node2; int BEtype;

                            fichier >> node1 >> node2;
                            if(DEBUG) std::cout << node1-one << " " << node2-one << std::endl;
    
                                 if(std::get<1>(gmsh_physicals[physTag])=="BCfile")             {BEtype=0;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="free")               {BEtype=1;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="wall")               {BEtype=2;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="outflow")            {BEtype=3;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedPressure")    {BEtype=4;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedVelocities")  {BEtype=5;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="axisymmetric_y")     {BEtype=6;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="axisymmetric_x")     {BEtype=7;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="BC_rec")                 {BEtype=10;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="free_rec")               {BEtype=11;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="wall_rec")               {BEtype=12;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="outflow_rec")            {BEtype=13;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedPressure_rec")    {BEtype=14;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedVelocities_rec")  {BEtype=15;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="axisymmetric_y_rec")     {BEtype=16;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="axisymmetric_x_rec")     {BEtype=17;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="recordingObject")        {BEtype=20;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="recObj")                 {BEtype=20;}
                            else   {std::cout << "WARNING: Boundary Element not defined!" << std::endl; BEtype=physTag;}

                            Line edge(id-one, Nodes[node1-one], Nodes[node2-one], BEtype, partitionTag-one);

                            if(DEBUG) edge.print();
                            if(BEtype==20) { // if it is an embedded recording object
                                Lines.push_back(edge);
                                CO_Lines_partitionTags.push_back(partitionTag-one);
                                m.CO_Lines++; // +1 to counter
                            } 
                            else {
                                BoundaryElements2d.push_back(edge);
                                KB_partitionTags.push_back(partitionTag-one);
                                m.KB++; // +1 to counter
                            }
                        }
                        if(elementType==2) // Triangle element
                        {   
                            // Nodes and boundaries type 
                            size_t node1, node2, node3; int Etype;

                            fichier >> node1 >> node2 >> node3;
                            if(DEBUG) std::cout << node1-one << " " << node2-one << " " << node3-one << std::endl;

                                 if(std::get<1>(gmsh_physicals[physTag])=="fluid")             {Etype=0;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid1")            {Etype=1;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid2")            {Etype=2;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid3")            {Etype=3;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid4")            {Etype=4;}
                            else   {std::cout << "WARNING: Fluid Element not defined!" << std::endl; Etype=physTag;}

                            Triangle element(id-one, Nodes[node1-one], Nodes[node2-one], Nodes[node3-one], Etype, partitionTag-one);

                            if(DEBUG) element.print();
                            Elements2d.push_back(element);
                            K_partitionTags.push_back(partitionTag-one);
                            m.K++; // +1 to counter
                        }
                        if(elementType==15) // Points (recording objects)
                        {   
                            // Single points are treated as virtual boundary elements
                            size_t node1; int BEtype;

                            fichier >> node1;
                            if(DEBUG) std::cout << node1-one << std::endl;

                                 if(std::get<1>(gmsh_physicals[physTag])=="recordingObject")   {BEtype=20;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="recObj")            {BEtype=20;}
                            else   {std::cout << "WARNING: Control Object not defined!" << std::endl; BEtype=physTag;}

                            Point point(id-one, Nodes[node1-one], BEtype, partitionTag-one);

                            if(DEBUG) point.print();
                            Points.push_back(point);
                            CO_Points_partitionTags.push_back(partitionTag-one);
                            m.CO_Points++; // +1 to counter
                        }
                    }
                }
                else if (phys_MAX_DIM==3) // We assume this is a 3-dimensional simulation problem.
                {
                    m.spatialDimension = 3;
                    m.Nfaces = 4;
                    if(DEBUG) std::cout << "Domain for 3-d test has been found." << std::endl;

                    // Read all types of elements and add them to their respective container
                    for (size_t i=0; i<num_elements; ++i)
                    {                        
                        fichier >> id >> elementType >> nTags;
                        if(DEBUG) std::cout << id-one << " " << elementType << " " << nTags << " ";

                        if(nTags==2)  // single domain
                        {  
                            fichier >> physTag >> geometryTag; num_PartitionTag=1; partitionTag=1;
                            if(DEBUG) std::cout << physTag << " " << geometryTag << " ";
                        } 
                        if(nTags==4) // partitioned domain
                        {   
                            fichier >> physTag >> geometryTag >> num_PartitionTag >> partitionTag;
                            if(DEBUG) std::cout << physTag << " " << geometryTag << " " << num_PartitionTag << " " << partitionTag-one << " ";
                        }
                        if(nTags>4) // partitioned domain
                        {   
                            fichier >> physTag >> geometryTag >> num_PartitionTag >> partitionTag;
                            if(DEBUG) std::cout << physTag << " " << geometryTag << " " << num_PartitionTag << " " << partitionTag-one << " ";

                            for (int i=0; i<(nTags-4); i++)
                            {
                                fichier >> negativeTag;
                                if(DEBUG) std::cout << negativeTag << " ";
                            }
                        }
                        if(elementType==1) // Lines (control objects)
                        {   
                            // Nodes and boundaries type 
                            size_t node1, node2; int BEtype;

                            fichier >> node1 >> node2;
                            if(DEBUG) std::cout << node1-one << " " << node2-one << std::endl;
    
                                 if(std::get<1>(gmsh_physicals[physTag])=="recordingObject")   {BEtype=20;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="recObj")            {BEtype=20;}
                            else   {std::cout << "WARNING: Boundary Element not defined!" << std::endl; BEtype=physTag;}

                            Line edge(id-one, Nodes[node1-one], Nodes[node2-one], BEtype, partitionTag-one);

                            if(DEBUG) edge.print();
                            Lines.push_back(edge);
                            CO_Lines_partitionTags.push_back(partitionTag-one);
                            m.CO_Lines++; // +1 to counter
                        }
                        if(elementType==2) // Triangle (Boundary) element
                        {   
                            // Nodes and boundaries type 
                            size_t node1, node2, node3; int BEtype;

                            fichier >> node1 >> node2 >> node3;
                            if(DEBUG) std::cout << node1-one << " " << node2-one << " " << node3-one << std::endl;
    
                                 if(std::get<1>(gmsh_physicals[physTag])=="BCfile")             {BEtype=0;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="free")               {BEtype=1;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="wall")               {BEtype=2;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="outflow")            {BEtype=3;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedPressure")    {BEtype=4;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedVelocities")  {BEtype=5;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="BC_gauss_spot")      {BEtype=6;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="BC_rec")                 {BEtype=10;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="free_rec")               {BEtype=11;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="wall_rec")               {BEtype=12;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="outflow_rec")            {BEtype=13;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedPressure_rec")    {BEtype=14;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="imposedVelocities_rec")  {BEtype=15;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="BC_gauss_spot_rec")      {BEtype=16;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="recordingObject")        {BEtype=20;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="recObj")                 {BEtype=20;}
                            else   {std::cout << "WARNING: Boundary Element not defined!" << std::endl; BEtype=physTag;}

                            Triangle face(id-one, Nodes[node1-one], Nodes[node2-one], Nodes[node3-one], BEtype, partitionTag-one);

                            if(DEBUG) face.print();
                            if(BEtype==20) { // if it is an embedded recording object
                                Faces.push_back(face);
                                CO_Faces_partitionTags.push_back(partitionTag-one);
                                m.CO_Faces++; // +1 to counter
                            } 
                            else {
                                BoundaryElements3d.push_back(face);
                                KB_partitionTags.push_back(partitionTag-one);
                                m.KB++; // +1 to counter
                            }
                        }
                        if(elementType==4) // Tetrahedron element
                        {   
                            // Nodes and boundaries type 
                            size_t node1, node2, node3, node4; int Etype;

                            fichier >> node1 >> node2 >> node3 >> node4;
                            if(DEBUG) std::cout << node1-one << " " << node2-one << " " << node3-one << " " << node4-one << std::endl;

                                 if(std::get<1>(gmsh_physicals[physTag])=="fluid")             {Etype=0;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid1")            {Etype=1;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid2")            {Etype=2;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid3")            {Etype=3;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="fluid4")            {Etype=4;}
                            else   {std::cout << "WARNING: Fluid Element not defined!" << std::endl; Etype=physTag;}

                            Tetrahedron element(id-one, Nodes[node1-one], Nodes[node2-one], Nodes[node3-one], Nodes[node4-one], Etype, partitionTag-one);

                            if(DEBUG) element.print();
                            Elements3d.push_back(element);
                            K_partitionTags.push_back(partitionTag-one);
                            m.K++; // +1 to counter
                        }
                        if(elementType==15) // Points (recording objects)
                        {   
                            // Single points are treated as virtual boundary elements
                            size_t node1; int BEtype;

                            fichier >> node1;
                            if(DEBUG) std::cout << node1-one << std::endl;

                                 if(std::get<1>(gmsh_physicals[physTag])=="recordingObject")   {BEtype=20;}
                            else if(std::get<1>(gmsh_physicals[physTag])=="recObj")            {BEtype=20;}
                            else   {std::cout << "WARNING: Control Object not defined!" << std::endl; BEtype=physTag;}

                            Point point(id-one, Nodes[node1-one], BEtype, partitionTag-one);

                            if(DEBUG) point.print();
                            Points.push_back(point);
                            CO_Points_partitionTags.push_back(partitionTag-one);
                            m.CO_Points++; // +1 to counter
                        }
                    }
                }
                else 
                {
                    // Other case is not available with the current solver
                    std::cout << " Error - The mesh is not consistent!" << std::endl; std::exit(-1);
                }
                if(DEBUG) std::cout << " boundary elements, KB: "<< m.KB << std::endl;
                if(DEBUG) std::cout << " domain elements, K : "<< m.K << std::endl;
                
                // count the number of elements/control objects in each partition
                std::map<size_t,size_t> Kmap; //auto& local_KP = KP;
                std::map<size_t,size_t> KBmap; //auto& local_KBP = KBP;
                std::map<size_t,size_t> CO_Points_map; 
                std::map<size_t,size_t> CO_Lines_map; 
                std::map<size_t,size_t> CO_Faces_map; 
                std::for_each( K_partitionTags.begin(), K_partitionTags.end(), [&Kmap]( size_t Tag ){ Kmap[Tag]++; } );
                std::for_each( KB_partitionTags.begin(), KB_partitionTags.end(), [&KBmap]( size_t Tag ){ KBmap[Tag]++; } );
                std::for_each( CO_Points_partitionTags.begin(), CO_Points_partitionTags.end(), [&CO_Points_map]( size_t Tag ){ CO_Points_map[Tag]++; } );
                std::for_each( CO_Lines_partitionTags.begin(), CO_Lines_partitionTags.end(), [&CO_Lines_map]( size_t Tag ){ CO_Lines_map[Tag]++; } );
                std::for_each( CO_Faces_partitionTags.begin(), CO_Faces_partitionTags.end(), [&CO_Faces_map]( size_t Tag ){ CO_Faces_map[Tag]++; } );

                // Total number of partitions in the k-elements domain
                m.num_partitions = Kmap.size();

                // Record the number of elements/control objects per partition
                for (size_t i=0; i<m.num_partitions; i++){
                    m.KP.push_back(Kmap[i]); // acumulate domain elements number
                    if (KBmap.count(i) > 0) { m.KBP.push_back(KBmap[i]);} // the domain has boundary data.
                    else {m.KBP.push_back(0);} // the domain DOESN'T has boundary data!
                    if (CO_Points_map.count(i) > 0) { m.CO_Points_P.push_back(CO_Points_map[i]);} // the domain has 3d point control objects.
                    else {m.CO_Points_P.push_back(0);} // the domain DOESN'T has 3d point control objects!
                    if (CO_Lines_map.count(i) > 0) { m.CO_Lines_P.push_back(CO_Lines_map[i]);} // the domain has 3d line control objects.
                    else {m.CO_Lines_P.push_back(0);} // the domain DOESN'T has 3d line control objects!
                    if (CO_Faces_map.count(i) > 0) { m.CO_Faces_P.push_back(CO_Faces_map[i]);} // the domain has 3d face control objects.
                    else {m.CO_Faces_P.push_back(0);} // the domain DOESN'T has 3d face control objects!
                };

                // Verify:
                if(DEBUG) std::cout << " mesh partitions : " << m.num_partitions << std::endl;
                if(DEBUG) std::cout << " partition | num domain elements" << std::endl;
                if(DEBUG) {for( auto p : Kmap ) std::cout << p.first << ' ' << p.second << std::endl;}
                if(DEBUG) std::cout << " partition | num boundary elements" << std::endl;
                if(DEBUG) {for( auto p : KBmap ) std::cout << p.first << ' ' << p.second << std::endl;}
                if(DEBUG) std::cout << " partition | num control points" << std::endl;
                if(DEBUG) {for( auto p : CO_Points_map ) std::cout << p.first << ' ' << p.second << std::endl;}
                if(DEBUG) std::cout << " partition | num control lines" << std::endl;
                if(DEBUG) {for( auto p : CO_Lines_map ) std::cout << p.first << ' ' << p.second << std::endl;}
                if(DEBUG) std::cout << " partition | num control faces" << std::endl;
                if(DEBUG) {for( auto p : CO_Faces_map ) std::cout << p.first << ' ' << p.second << std::endl;}

                // finish reading the last line
                std::getline(fichier, ligne);

                // Verify that is the actual end of the section
                std::getline(fichier, ligne);
                if(ligne.find("$EndElements") != static_cast<std::string::size_type>(0))
                {
                    std::cout << " Error - Something is wrong with the elements!" << std::endl; std::exit(-1);
                } 
            }
            // If we are at the end of the file 
            else if (fichier.eof())
            {
                break;
            }
        }
    } else {
        std::cout << "Could not open file: " << mesh_file << std::endl; std::exit(-1);
    }  
    
    if (phys_MAX_DIM==2)
    {
        // Assume triangle elements for the domain
        size_t Nfaces = 3; // number of faces
        size_t D = 2; // physical dimension
        size_t n; // object counter

        //--------------------------------------------
        //  3 - BUILD VERTEX (V) ARRAY
        //--------------------------------------------

        // contains the vertex the Skeleton mesh nodes!
        std::vector<size_t>dims_V={ Nodes.size(),D}; m.V.allocate(dims_V);

        n=0;
        for ( const auto node : Nodes ) {
            m.V(n,0)=node.get_x();
            m.V(n,1)=node.get_y();
            n++;
        }

        //-------------------------------------------------
        //  4 - ELEMENT TO VERTEX COORDINATES (EToVc) MAP
        //-------------------------------------------------

        // contains the vextex coordinates of each triangle element 
        std::vector<size_t>dims_EToVc={ Elements2d.size(),D*Nfaces}; m.EToVc.allocate(dims_EToVc);

        n=0;
        for ( const auto element : Elements2d ) {
            // Node 1
            m.EToVc(n,0)=element.get_N1x();
            m.EToVc(n,1)=element.get_N1y();
            // Node 2
            m.EToVc(n,2)=element.get_N2x();
            m.EToVc(n,3)=element.get_N2y();
            // Node 3
            m.EToVc(n,4)=element.get_N3x();
            m.EToVc(n,5)=element.get_N3y();
            // Save element type
            m.Etype.push_back(element.get_ElementType());
            n++;
        }

        //--------------------------------------------
        //  5 - ELEMENT TO VERTEX (EToV) MAP
        //--------------------------------------------

        // Element to vertex map for the triangles element
        std::vector<size_t>dims_EToV={ Elements2d.size(),Nfaces}; m.EToV.allocate(dims_EToV);

        n=0;
        for ( const auto element : Elements2d ) {
            // Element To Vertex (EToV) connectivity
            m.EToV(n,0)=element.get_N1Id();
            m.EToV(n,1)=element.get_N2Id();
            m.EToV(n,2)=element.get_N3Id();
            n++;
        }

        //------------------------------------------------
        //  6 - BOUNDARY ELEMENTS TO VERTEX (BEToV) MAP
        //------------------------------------------------

        // contains the vertex, the type of the faces at the boundary
        std::vector<size_t>dims_BEToV={ BoundaryElements2d.size(),Nfaces}; m.BEToV.allocate(dims_BEToV);

        n=0;
        for ( const auto edge : BoundaryElements2d ) {
            // Boundary Element To Vertex (BEToV) connectivity
            m.BEToV(n,0) = edge.get_N1Id();
            m.BEToV(n,1) = edge.get_N2Id();
            m.BEToV(n,2) = edge.get_EdgeType();
            n++;
        }

        //------------------------------------------------
        //  7 - ELEMENTS TO PARTITION (EToP) MAP
        //------------------------------------------------

        // contains the vertex, the type of the faces at the boundary
        std::vector<size_t>dims_EToP={ Elements2d.size(),2}; m.EToP.allocate(dims_EToP);

        n=0;
        for ( const auto element : Elements2d ) {
            // Element To Partition (EToP) map
            m.EToP(n,0) = element.get_Partition();
            m.EToP(n,1) = n; // Element Index 
            n++;
        }

        //---------------------------------------------------
        //  8 - BOUNDARY ELEMENTS TO PARTITION (BEToP) MAP
        //---------------------------------------------------

        // contains the vertex, the type of the faces at the boundary
        std::vector<size_t>dims_BEToP={ BoundaryElements2d.size(),2}; m.BEToP.allocate(dims_BEToP);

        n=0;
        for ( const auto edge : BoundaryElements2d ) {
            // Element To Partition (EToP) map
            m.BEToP(n,0) = edge.get_Partition();
            m.BEToP(n,1) = n; // Element Index 
            n++;
        }

        //---------------------------------------------------
        //  9 - BUILD VERTEX OF CONTROL OBJECTS
        //---------------------------------------------------

        // contains the vertex of Points
        std::vector<size_t>dims_PointsToV={ Points.size(),1+1}; m.PointsToV.allocate(dims_PointsToV);

        n=0;
        for ( const auto point : Points ) {
            m.PointsToV(n,0)=point.get_Partition(); // Partition
            m.PointsToV(n,1)=point.get_N1Id(); // Node 1
            n++;
        }

        // contains the vertex of Lines
        std::vector<size_t>dims_LinesToV={ Lines.size(),D+1}; m.LinesToV.allocate(dims_LinesToV);

        n=0;
        for ( const auto line : Lines ) {
            m.LinesToV(n,0)=line.get_Partition(); // Partition
            m.LinesToV(n,1)=line.get_N1Id(); // Node 1
            m.LinesToV(n,2)=line.get_N2Id(); // Node 2
            n++;
        }

    }
    if (phys_MAX_DIM==3)
    {
        // Assume tetrahedron elements for the domain
        size_t Nfaces = 4; // number of faces
        size_t D = 3; // physical dimension
        size_t n; // object counter

        //--------------------------------------------
        //  3 - BUILD VERTEX (V) ARRAY
        //--------------------------------------------

        // contains the vertex the Skeleton mesh nodes!
        std::vector<size_t>dims_V={ Nodes.size(),D}; m.V.allocate(dims_V);

        n=0;
        for ( const auto node : Nodes ) {
            m.V(n,0)=node.get_x();
            m.V(n,1)=node.get_y();
            m.V(n,2)=node.get_z();
            n++;
        }

        //-------------------------------------------------
        //  4 - ELEMENT TO VERTEX COORDINATES (EToVc) MAP
        //-------------------------------------------------

        // contains the vextex coordinates of each triangle element 
        std::vector<size_t>dims_EToVc={ Elements3d.size(),D*Nfaces}; m.EToVc.allocate(dims_EToVc);

        n=0;
        for ( const auto element : Elements3d ) {
            // Node 1
            m.EToVc(n,0)=element.get_N1x();
            m.EToVc(n,1)=element.get_N1y();
            m.EToVc(n,2)=element.get_N1z();
            // Node 2
            m.EToVc(n,3)=element.get_N2x();
            m.EToVc(n,4)=element.get_N2y();
            m.EToVc(n,5)=element.get_N2z();
            // Node 3
            m.EToVc(n,6)=element.get_N3x();
            m.EToVc(n,7)=element.get_N3y();
            m.EToVc(n,8)=element.get_N3z();
            // Node 4
            m.EToVc(n,9)=element.get_N4x();
            m.EToVc(n,10)=element.get_N4y();
            m.EToVc(n,11)=element.get_N4z();
            // Save element type
            m.Etype.push_back(element.get_ElementType());
            n++;
        }

        //--------------------------------------------
        //  5 - ELEMENT TO VERTEX (EToV) MAP
        //--------------------------------------------

        // Element to vertex map for the triangles element
        std::vector<size_t>dims_EToV={ Elements3d.size(), Nfaces}; m.EToV.allocate(dims_EToV);

        n=0;
        for ( const auto element : Elements3d ) {
            // Element To Vertex (EToV) connectivity
            m.EToV(n,0)=element.get_N1Id();
            m.EToV(n,1)=element.get_N2Id();
            m.EToV(n,2)=element.get_N3Id();
            m.EToV(n,3)=element.get_N4Id();
            n++;
        }

        //------------------------------------------------
        //  6 - BOUNDARY ELEMENTS TO VERTEX (BEToV) MAP
        //------------------------------------------------

        // contains the vertex, the type of the faces at the boundary
        std::vector<size_t>dims_BEToV={ BoundaryElements3d.size(),Nfaces}; m.BEToV.allocate(dims_BEToV);

        n=0;
        for ( const auto face : BoundaryElements3d ) {
            // Boundary Element To Vertex (BEToV) connectivity
            m.BEToV(n,0) = face.get_N1Id();
            m.BEToV(n,1) = face.get_N2Id();
            m.BEToV(n,2) = face.get_N3Id();
            m.BEToV(n,3) = face.get_FaceType();
            n++;
        }

        //------------------------------------------------
        //  7 - ELEMENTS TO PARTITION (EToP) MAP
        //------------------------------------------------

        // contains the vertex, the type of the faces at the boundary
        std::vector<size_t>dims_EToP={ Elements3d.size(),2}; m.EToP.allocate(dims_EToP);

        n=0;
        for ( const auto element : Elements3d ) {
            // Element To Partition (EToP) map
            m.EToP(n,0) = element.get_Partition();
            m.EToP(n,1) = n; // Element Index 
            n++;
        }

        //---------------------------------------------------
        //  8 - BOUNDARY ELEMENTS TO PARTITION (BEToP) MAP
        //---------------------------------------------------

        // contains the vertex, the type of the faces at the boundary
        std::vector<size_t>dims_BEToP={ BoundaryElements3d.size(),2}; m.BEToP.allocate(dims_BEToP);

        n=0;
        for ( const auto face : BoundaryElements3d ) {
            // Element To Partition (EToP) map
            m.BEToP(n,0) = face.get_Partition();
            m.BEToP(n,1) = n; // Element Index 
            n++;
        }

        //---------------------------------------------------
        //  9 - BUILD VERTEX OF CONTROL OBJECTS
        //---------------------------------------------------

        // contains the vertex of Points
        std::vector<size_t>dims_PointsToV={ Points.size(),1+1}; m.PointsToV.allocate(dims_PointsToV);

        n=0;
        for ( const auto point : Points ) {
            m.PointsToV(n,0)=point.get_Partition(); // Partition
            m.PointsToV(n,1)=point.get_N1Id(); // Node 1
            n++;
        }

        // contains the vertex of Lines
        std::vector<size_t>dims_LinesToV={ Lines.size(),2+1}; m.LinesToV.allocate(dims_LinesToV);

        n=0;
        for ( const auto line : Lines ) {
            m.LinesToV(n,0)=line.get_Partition(); // Partition
            m.LinesToV(n,1)=line.get_N1Id(); // Node 1
            m.LinesToV(n,2)=line.get_N2Id(); // Node 2
            n++;
        }

        // contains the vertex of triangular faces
        std::vector<size_t>dims_FacesToV={ Faces.size(),D+1}; m.FacesToV.allocate(dims_FacesToV);

        n=0;
        for ( const auto face : Faces ) {
            m.FacesToV(n,0)=face.get_Partition(); // Partition
            m.FacesToV(n,1)=face.get_N1Id(); // Node 1            
            m.FacesToV(n,2)=face.get_N2Id(); // Node 2            
            m.FacesToV(n,3)=face.get_N3Id(); // Node 3
            n++;
        }
    }
    // Lastly, if everything is normal ...
    return m;
}

#endif