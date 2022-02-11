function mesh = GMSHv4(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Extract entities contained in a single GMSH file in format v4.1 
%
%      Coded by Manuel A. Diaz @ Univ-Poitiers | Pprime, 2022.01.21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example call: mesh = GMSHv4('filename.msh');
%
% Output: 
%   mesh.Dim:  mesh dimensionality, D=2,3.
%   mesh.V :   vertices/nodes coordinates.
%   mesh.E1:   2-node line.
%   mesh.E2:   3-node triangle.
%   mesh.E3:   4-node quadrangle.
%   mesh.E4:   4-node tetrahedron.
%   mesh.E5:   8-node hexahedron.
%   mesh.E6:   6-node prism.
%   mesh.E7:   5-node pyramid.
%   mesh.E8:   3-node second order line.
%   mesh.E9:   6-node second order triangle.
%   mesh.E10:  9-node second order quadrangle.
%   mesh.E11: 10-node second order tetrahedron.
%   mesh.E12: 27-node second order hexahedron.
%   mesh.E13: 18-node second order prism.
%   mesh.E14: 14-node second order pyramid.
%   mesh.E15:  1-node point.
%   mesh.E16:  8-node second order quadrangle.
%   mesh.E17: 20-node second order hexahedron.
%   mesh.E18: 15-node second order prism.
%   mesh.E19: 13-node second order pyramid.
%   mesh.E20:  9-node third order incomplete triangle.
%   mesh.E21: 10-node third order triangle.
%   mesh.E22: 12-node fourth order incomplete triangle.
%   mesh.E23: 15-node fourth order triangle.
%   mesh.E24: 15-node fifth order incomplete triangle.
%   mesh.E25: 21-node fifth order complete triangle.
%   mesh.E26:  4-node third order edge.
%   mesh.E27:  5-node fourth order edge.
%   mesh.E28:  6-node fifth order edge.
%   mesh.E29: 20-node third order tetrahedron.
%   mesh.E30: 35-node fourth order tetrahedron.
%   mesh.E31: 56-node fifth order tetrahedron.
%   mesh.E92: 64-node third order hexahedron.
%   mesh.E93:125-node fourth order hexahedron.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read all sections
DEBUG = true;

file = fileread(filename);
% Extract strings between:
MeshFormat    = extractBetween(file,'$MeshFormat','$EndMeshFormat');
PhysicalNames = extractBetween(file,'$PhysicalNames','$EndPhysicalNames');
Entities      = extractBetween(file,'$Entities','$EndEntities');
PartEntities  = extractBetween(file,'$PartitionedEntities','$EndPartitionedEntities');
Nodes         = extractBetween(file,'$Nodes','$EndNodes');
Elements      = extractBetween(file,'$Elements','$EndElements');

% Split data lines into cells
cells_MF   = splitlines(MeshFormat);
cells_PN   = splitlines(PhysicalNames);
cells_Ent  = splitlines(Entities);
cells_PEnt = splitlines(PartEntities);
cells_N    = splitlines(Nodes);
cells_E    = splitlines(Elements);

% Delete emptly cells 
cells_MF   = cells_MF(not(cellfun('isempty',cells_MF)));
cells_PN   = cells_PN(not(cellfun('isempty',cells_PN)));
cells_Ent  = cells_Ent(not(cellfun('isempty',cells_Ent)));
cells_PEnt = cells_PEnt(not(cellfun('isempty',cells_PEnt)));
cells_N    = cells_N(not(cellfun('isempty',cells_N)));
cells_E    = cells_E(not(cellfun('isempty',cells_E)));

% Sanity check
if cellfun('isempty',cells_MF ), error('Error - Wrong File Format!'); end
if cellfun('isempty',cells_PN ), error('Error - No Physical names!'); end
if cellfun('isempty',cells_Ent), error('Error - No Entities found!'); end
if cellfun('isempty',cells_N  ), error('Error - Nodes are missing!'); end
if cellfun('isempty',cells_E  ), error('Error - No elements found!'); end

%% Identify critical data within each section:

%********************%
% 1. Read MeshFormat
%********************%
line_data = sscanf(cells_MF{1},'%f %d %d');
mesh.version   = line_data(1);	% 4.1 is expected
mesh.file_type = line_data(2);	% 0:ASCII or 1:Binary
mesh.mode      = line_data(3);	% 1 in binary mode to detect endianness 

% Sanity check
if (mesh.version ~= 4.1), error('Error - Expected mesh format v4.1'); end
if (mesh.file_type ~= 0), error('Error - Binary file not allowed'); end

%***********************%
% 2. Read PhysicalNames
%***********************%
numPhysicalNames = sscanf(cells_PN{1},'%d');
for n = 1:numPhysicalNames
   line_data = sscanf(cells_PN{n+1},'%d %d');
   mesh.physicalNames(n).dimension = line_data(1);
   mesh.physicalNames(n).tag = line_data(2);
   line_data = regexp(cells_PN{n+1},'(?<=")[^"]+(?=")','match');
   mesh.physicalNames(n).name = line_data{1};
end

%***********************%
% 3. Read Entities
%***********************%
l=1; % line counter
line_data = sscanf(cells_Ent{l},'%d %d %d %d');
numPoints = line_data(1);
numCurves = line_data(2);
numSurfaces= line_data(3);
numVolumes = line_data(4);
l=2; % line counter
%mesh.entities.points
for i = 1:numPoints
    [pointTag,pointPhysicalTags] = get_entity(cells_Ent{l},'node');
    disp([pointTag,pointPhysicalTags]);
    l = l+1; % update line counter
end
%mesh.entities.curves
for i = 1:numCurves
    [curveTag,curvePhysicalTags] = get_entity(cells_Ent{l},'curve');
    disp([curveTag,curvePhysicalTags]);
    l = l+1; % update line counter
end
%mesh.entities.surfaces
for i = 1:numSurfaces
    [surfaceTag,surfacePhysicalTags] = get_entity(cells_Ent{l},'surface');
    disp([surfaceTag,surfacePhysicalTags]);
    l = l+1; % update line counter
end
%mesh.entities.volumes
for i = 1:numVolumes
    [volumeTag,volumePhysicalTags] = get_entity(cells_Ent{l},'volume');
    disp([volumeTag,volumePhysicalTags]);
    l = l+1; % update line counter
end

%******************************%
% 4. Read Partitioned Entities
%******************************%
l=1; numPartitions = sscanf(cells_PEnt{l},'%d');
l=2; numGhostEntities = sscanf(cells_PEnt{l},'%d'); % not important for the moment!
l=3; line_data = sscanf(cells_PEnt{l},'%d');
numPoints  = line_data(1);
numCurves  = line_data(2);
numSurfaces= line_data(3);
numVolumes = line_data(4);
l=4; % line counter
%mesh.partitionedEntities.points
for i = 1:numPoints
    [pointTag,pointPartTags,pointPhysicalTags] = get_partitionedEntity(cells_PEnt{l},'node');
    disp([pointTag,pointPartTags,pointPhysicalTags]);
    l = l+1; % update line counter
end
%mesh.partitionedEntities.curves
for i = 1:numCurves
    [curveTag,curvePartTags,curvePhysicalTags] = get_partitionedEntity(cells_PEnt{l},'curve');
    disp([curveTag,curvePartTags,curvePhysicalTags]);
    l = l+1; % update line counter
end
%mesh.partitionedEntities.surfaces
for i = 1:numSurfaces
    [surfaceTag,surfacePartTags,surfacePhysicalTags] = get_partitionedEntity(cells_PEnt{l},'surface');
    disp([surfaceTag,surfacePartTags,surfacePhysicalTags]);
    l = l+1; % update line counter
end
%mesh.partitionedEntities.volumes
for i = 1:numVolumes
    [volumeTag,volumePartTags,volumePhysicalTags] = get_partitionedEntity(cells_PEnt{l},'volume');
    disp([volumeTag,volumePartTags,volumePhysicalTags]);
    l = l+1; % update line counter
end

%***********************%
% Read Nodes
%***********************%
l=1; % line counter
line_data = sscanf(cells_N{l},'%d %d %d %d');
numEntityBlocks = line_data(1);
numNodes        = line_data(2);
minNodeTag      = line_data(3);
maxNodeTag      = line_data(4);

% allocate space for nodal data
mesh.V = zeros(numNodes,3); % [x,y,z] 

% Note: entityBlock can be read in parallel !
for ent = 1:numEntityBlocks
    l = l+1; % update line counter
    line_data = sscanf(cells_N{l},'%d %d %d %d');
    entityDim = line_data(1);
    entityTag = line_data(2);
    %parametric = line_data(3); % not used for the moment.
    numNodesInBlock = line_data(4);
    %
    nodeTag = zeros(1,numNodesInBlock); % nodeTag
    for i=1:numNodesInBlock
        l = l+1; % update line counter
        nodeTag(i) = sscanf(cells_N{l},'%d');
    end
    %
    for i=1:numNodesInBlock
        l = l+1; % update line counter
        mesh.V(nodeTag(i),:) = sscanf(cells_N{l},'%f %f %f'); % [x(i),y(i),z(i)]
    end
end

if DEBUG 
    x = double(mesh.V(:,1));
    y = double(mesh.V(:,2));
    z = double(mesh.V(:,3));
    scatter3(x,y,z,'.k'); 
    hold on
    for i=1:numNodes
        text(x(i),y(i),z(i),num2str(i));
    end
    hold off
end

%***********************%
% Read Elements
%***********************%
l=1; % line counter
line_data = sscanf(cells_E{l},'%d %d %d %d');
numEntityBlocks = line_data(1);
numElements     = line_data(2);
minElementsTag  = line_data(3);
maxElementsTag  = line_data(4);

% Allocate space for Elements data
mesh.E1 =sparse(numElements,2);
mesh.E2 =sparse(numElements,3);
mesh.E4 =sparse(numElements,4);
mesh.E15=sparse(numElements,1);

for ent = 1:numEntityBlocks
    l = l+1; % update line counter
    line_data = sscanf(cells_E{l},'%d %d %d %d');
    entityDim = line_data(1); % entity spatial dimension
    entityTag = line_data(2); % Subgroup Tag number
    elementType = line_data(3);
    numElementsInBlock = line_data(4);
    %
    elementTag = zeros(1,numElementsInBlock); % elementTag
    for i=1:numElementsInBlock
        l = l+1; % update line counter
        line_data = sscanf(cells_E{l},'%d %d %d %d');
        elementTag(i) = line_data(1);
        switch elementType
            case 1, mesh.E1 (elementTag(i),:) = line_data(2:end);
            case 2, mesh.E2 (elementTag(i),:) = line_data(2:end);
            case 4, mesh.E4 (elementTag(i),:) = line_data(2:end);
            case 15,mesh.E15(elementTag(i),:) = line_data(2);
            otherwise, error('element not in list');
        end
    end
end
%
end % GMSHv4 read function

% Get single entity information:
function [entityTag,physicalTags] = get_entity(str_line,type)
    
    vector = str2double(regexp(str_line,'-?[\d.]+(?:e-?\d+)?','match'));

    switch type
        case 'node'
            % 1. get entityTag
            entityTag = vector(1);

            % 2. get entity coordinates % not needed
            %minX(i) = line_data(2); 
            %minY(i) = line_data(3);
            %minZ(i) = line_data(4);

            % 3. get physical tag associated
            numPhysicalTags = vector(5);
            physicalTags = zeros(1,numPhysicalTags);
            for j=1:numPhysicalTags
                physicalTags(j) = vector(5+j);
            end

        otherwise
            % 1. get entityTag
            entityTag = vector(1);
        
            % 2. get entity limits (for visualization) % not needed
            %minX(i) = line_data(2); 
            %minY(i) = line_data(3);
            %minZ(i) = line_data(4);
            %maxX(i) = line_data(5);
            %maxY(i) = line_data(6);
            %maxZ(i) = line_data(7);
            
            % 3. get physical tag associated
            numPhysicalTags = vector(8);
            physicalTags = zeros(1,numPhysicalTags);
            for j=1:numPhysicalTags
                physicalTags(j) = vector(8+j);
            end
        
            % 4. get tags of subentities that define it. % not needed
            %numBoudingEntities = line_data(9+j);
            %entitiesTags = zeros(1,numBoudingEntities);
            %for k=1:numBoudingEntities
            %   entitiesTags(k) = line_data(9+j+k);
            %end
    end
end

% Get single partitioned entity information:
function [entityTag,partitionTags,physicalTags] = get_partitionedEntity(str_line,type)
    
    vector = str2double(regexp(str_line,'-?[\d.]+(?:e-?\d+)?','match'));

    switch type
        case 'node'
            % 1. get entityTag
            entityTag = vector(1);

            % 2. get parent dimention and tag % no needed
            %parentDim(i) = vector(2);
            %parentTag(i) = vector(3);
            
            numPartitionTags = vector(4);
            partitionTags = zeros(1,numPartitionTags);
            for j=1:numPartitionTags
                partitionTags(j) = vector(4+j);
            end

            % 3. get entity coordinates % not needed
            %minX(i) = line_data(5+j); 
            %minY(i) = line_data(6+j);
            %minZ(i) = line_data(7+j);

            % 4. get physical tag associated
            numPhysicalTags = vector(8+j);
            physicalTags = zeros(1,numPhysicalTags);
            for k=1:numPhysicalTags
                physicalTags(k) = vector(8+j+k);
            end

        otherwise
            % 1. get entityTag
            entityTag = vector(1);
        
            % 2. get parent dimention and tag % no needed
            %parentDim(i) = vector(2);
            %parentTag(i) = vector(3);
            
            numPartitionTags = vector(4);
            partitionTags = zeros(1,numPartitionTags);
            for j=1:numPartitionTags
                partitionTags(j) = vector(4+j);
            end

            % 3. get entity coordinates % not needed
            %minX(i) = line_data( 5+j); 
            %minY(i) = line_data( 6+j);
            %minZ(i) = line_data( 7+j);
            %maxX(i) = line_data( 8+j); 
            %maxY(i) = line_data( 9+j);
            %maxZ(i) = line_data(10+j);
            
            % 4. get physical tag associated
            numPhysicalTags = vector(11+j);
            physicalTags = zeros(1,numPhysicalTags);
            for k=1:numPhysicalTags
                physicalTags(k) = vector(11+j+k);
            end
        
            % 5. get tags of subentities that define it. % not needed
            %numBoudingEntities = line_data(12+j+k);
            %entitiesTags = zeros(1,numBoudingEntities);
            %for l=1:numBoudingEntities
            %   entitiesTags(l) = line_data(12+j+k+l);
            %end
    end
end