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

f = fileread(filename);

% Extract strings between:
MeshFormat           = extractBetween(f,'$MeshFormat','$EndMeshFormat');
PhysicalNames        = extractBetween(f,'$PhysicalNames','$EndPhysicalNames');
Entities             = extractBetween(f,'$Entities','$EndEntities');
PartitionedEntities  = extractBetween(f,'$PartitionedEntities','$EndPartitionedEntities');
Nodes                = extractBetween(f,'$Nodes','$EndNodes');
Elements             = extractBetween(f,'$Elements','$EndElements');
Periodic             = extractBetween(f,'$Periodic','$EndPeriodic');
GhostElements        = extractBetween(f,'$GhostElements','$EndGhostElements');
Parametrizations     = extractBetween(f,'$GhostParametrizations','$EndGhostParametrizations');
NodeData             = extractBetween(f,'$NodeData','$EndNodeData');
ElementData          = extractBetween(f,'$ElementData','$EndElementData');
InterpolationMatrices= extractBetween(f,'$InterpolationMatrices','$EndInterpolationMatrices');

% Split data lines into cells
cells_MF    = splitlines(MeshFormat);
cells_PN    = splitlines(PhysicalNames);
cells_Ent   = splitlines(Entities);
cells_ParEnt= splitlines(PartitionedEntities);
cells_N     = splitlines(Nodes);
cells_E     = splitlines(Elements);
cells_P     = splitlines(Periodic);
cells_GE    = splitlines(GhostElements);
cells_Param = splitlines(Parametrizations);
cells_ND    = splitlines(NodeData);
cells_ED    = splitlines(ElementData);
cells_IM    = splitlines(InterpolationMatrices);

% Delete emptly cells
cells_MF    = cells_MF(not(cellfun('isempty',cells_MF)));
cells_PN    = cells_PN(not(cellfun('isempty',cells_PN)));
cells_Ent   = cells_Ent(not(cellfun('isempty',cells_Ent)));
cells_ParEnt= cells_ParEnt(not(cellfun('isempty',cells_ParEnt)));
cells_N     = cells_N(not(cellfun('isempty',cells_N)));
cells_E     = cells_E(not(cellfun('isempty',cells_E)));
cells_P     = cells_P(not(cellfun('isempty',cells_P)));
cells_GE    = cells_GE(not(cellfun('isempty',cells_GE)));
cells_Param = cells_Param(not(cellfun('isempty',cells_Param)));
cells_ND    = cells_ND(not(cellfun('isempty',cells_ND)));
cells_ED    = cells_ED(not(cellfun('isempty',cells_ED)));
cells_IM    = cells_IM(not(cellfun('isempty',cells_IM)));

%% Identify critical data within each section:

% 1. Get MeshFormat
line_data = sscanf(cells_MF{1},'%f %d %d');
mesh.version   = line_data(1);	% 4.1 is expected
mesh.file_type = line_data(2);	% 0:ASCII or 1:Binary
mesh.mode      = line_data(3);	% 1 in binary mode to detect endianness 

% 2. Get PhysicalNames
numPhysicalNames = sscanf(cells_PN{1},'%d');
for n = 1:numPhysicalNames
   line_data = sscanf(cells_PN{n+1},'%d %d');
   mesh.physicalNames(n).dimension = line_data(1);
   mesh.physicalNames(n).tag = line_data(2);
   line_data = regexp(cells_PN{n+1},'(?<=")[^"]+(?=")','match');
   mesh.physicalNames(n).name = line_data{1};
end

% 3. Get Entities
l=1; % line counter
line_data = sscanf(cells_Ent{l},'%d %d %d %d');
numPoints = line_data(1);
numCurves = line_data(2);
numSurfaces= line_data(3);
numVolumes = line_data(4);
l=2; % line counter
%mesh.entities.points
for i = 1:numPoints
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    pointTag(i) = line_data(1);
    X(i) = line_data(2);
    Y(i) = line_data(3);
    Z(i) = line_data(4);
    numPhysicalTags = line_data(5);
    for j=1:numPhysicalTags
        physicalTag(i) = line_data(5+j);
    end
    l = l+1; % update line counter
end
%mesh.entities.curves
for i = 1:numCurves
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    curveTag(i) = line_data(1);
    minX(i) = line_data(2);
    minY(i) = line_data(3);
    minZ(i) = line_data(4);
    maxX(i) = line_data(5);
    maxY(i) = line_data(6);
    maxZ(i) = line_data(7);
    numPhysicalTags = line_data(8);
    for j=1:numPhysicalTags
        physicalTags(i) = line_data(8+j);
    end
    numBoundingPoints = line_data(9+j);
    for k=1:numBoundingPoints
        pointTags(i,k) = line_data(9+j+k);
    end
    l = l+1; % update line counter
end
%mesh.entities.surfaces
for i = 1:numSurfaces
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    surfaceTag(i) = line_data(1);
    minX(i) = line_data(2);
    minY(i) = line_data(3);
    minZ(i) = line_data(4);
    maxX(i) = line_data(5);
    maxY(i) = line_data(6);
    maxZ(i) = line_data(7);
    numPhysicalTags = line_data(8);
    for j=1:numPhysicalTags
        physicalTag(i) = line_data(9+j);
    end
    numBoundingCurves = line_data(10);
    for j=1:numBoundingCurves
        curveTag(i) = line_data(11);
    end
    l = l+1; % update line counter
end
%mesh.entities.volumes
for i = 1:numVolumes
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    volumeTag(i) = line_data(1);
    minX(i) = line_data(2);
    minY(i) = line_data(3);
    minZ(i) = line_data(4);
    maxX(i) = line_data(5);
    maxY(i) = line_data(6);
    maxZ(i) = line_data(7);
    numPhysicalTags = line_data(8);
    for j=1:numPhysicalTags
        physicalTag = line_data(9);
    end
    numBoundingSurfaces = line_data(10);
    for j=1:numBoundingSurfaces
        surfaceTag = line_data(11);
    end
    l = l+1; % update line counter
end

% 4. Get Partitioned Entities
l=1; numPartitions = sscanf(cells_ParEnt{l},'%d');
l=2; numGhostEntities = sscanf(cells_ParEnt{l},'%d'); % not important for the moment!
l=3; line_data = sscanf(cells_ParEnt{l},'%d');
numPoints  = line_data(1);
numCurves  = line_data(2);
numSurfaces= line_data(3);
numVolumes = line_data(4);
l=4; % line counter
%mesh.partitionedEntities.points
for i = 1:numPoints
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    pointTag(i) = line_data(1);
    parentDim(i) = line_data(2);
    parentTag(i) = line_data(3);
    numPartitions(i) = line_data(4);
    for j=1:numPartitions(i)
        partitionTag(i) = line_data(4+j);
    end
    %coords(i) = line_data(2:4);
    numPhysicalTags(i) = line_data(4+j+1);
    for j=1:numPhysicalTags(i)
        physicalTags(i) = line_data(6);
    end
    l = l+1; % update line counter
end
%mesh.partitionedEntities.curves
for i = 1:numCurves
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    l = l+1; % update line counter
end
%mesh.partitionedEntities.surfaces
for i = 1:numSurfaces
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    l = l+1; % update line counter
end
%mesh.partitionedEntities.volumes
for i = 1:numVolumes
    line_data = str2double(regexp(cells_Ent{l},'-?\d+\.?\d*','match'));
    l = l+1; % update line counter
end

% Get Nodes
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
    parametric = line_data(3); % we are not interested parametric entities for the moment.
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

% Get Elements
l=1; % line counter
line_data = sscanf(cells_E{l},'%d %d %d %d');
numEntityBlocks = line_data(1);
numElements     = line_data(2);
minElementsTag  = line_data(3);
maxElementsTag  = line_data(4);

% Allocate space for Elements data
mesh.E1 =sparse(numElements,2);
mesh.E2 =sparse(numElements,3);
mesh.E15=sparse(numElements,1);

for ent = 1:numEntityBlocks
    l = l+1; % update line counter
    line_data = sscanf(cells_E{l},'%d %d %d %d');
    entityDim = line_data(1); % entity spatial dimension
    entityTag = line_data(2); disp(entityTag); % Subgroup Tag number
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
            case 15,mesh.E15(elementTag(i),:) = line_data(2);
            otherwise, error('element not in list');
        end
    end
end
%
end % GMSHv4 read function