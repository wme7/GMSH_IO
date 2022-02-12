function [V,VE,SE,LE,PE,mapPhysNames,info] = GMSHparserV4(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Extract entities contained in a single GMSH file in format v4.1 
%
%      Coded by Manuel A. Diaz @ Pprime | Univ-Poitiers, 2022.01.21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example call: GMSHparserV4('filename.msh')
%
% Output:
%     V: the vertices (nodes coordinates) -- (Nx3) array
%    VE: volumetric elements (tetrahedrons) -- structure
%    SE: surface elements (triangles,quads) -- structure
%    LE: curvilinear elements (lines/edges) -- structure
%    PE: point elements (singular vertices) -- structure
%    mapPhysNames: maps phys.tag --> phys.name  -- map structure
%    info: version, format, endian test -- structure
%
% Note: This parser is designed to capture the following elements:
%
% Point:                Line:                   Triangle:
%
%        v                                              v
%        ^                                              ^
%        |                       v                      |
%        +----> u                ^                      2
%       0                        |                      |`\
%                                |                      |  `\
%                          0-----+-----1 --> u          |    `\
%                                                       |      `\
% Tetrahedron:                                          |        `\
%                                                       0----------1 --> u
%                    v
%                   ,
%                  /
%               2
%             ,/|`\                    Based on the GMSH guide 4.9.4
%           ,/  |  `\                  This are lower-order elements 
%         ,/    '.   `\                identified as:
%       ,/       |     `\                  E-1 : 2-node Line 
%     ,/         |       `\                E-2 : 3-node Triangle
%    0-----------'.--------1 --> u         E-4 : 4-node tetrahedron
%     `\.         |      ,/                E-15: 1-node point
%        `\.      |    ,/
%           `\.   '. ,/                Other elements can be added to 
%              `\. |/                  this parser by modifying the 
%                 `3                   Read Elements stage.
%                    `\.
%                       ` w            Happy coding ! M.D. 02/2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read all sections

file = fileread(filename);
% Extract strings between:
MeshFormat    = extractBetween(file,['$MeshFormat',newline],[newline,'$EndMeshFormat']);
PhysicalNames = extractBetween(file,['$PhysicalNames',newline],[newline,'$EndPhysicalNames']);
Entities      = extractBetween(file,['$Entities',newline],[newline,'$EndEntities']);
PartEntities  = extractBetween(file,['$PartitionedEntities',newline],[newline,'$EndPartitionedEntities']);
Nodes         = extractBetween(file,['$Nodes',newline],[newline,'$EndNodes']);
Elements      = extractBetween(file,['$Elements',newline],[newline,'$EndElements']);

% Sanity check
if isempty(MeshFormat),    error('Error - Wrong File Format!'); end
if isempty(PhysicalNames), error('Error - No Physical names!'); end
if isempty(Entities),      error('Error - No Entities found!'); end
if isempty(Nodes),         error('Error - Nodes are missing!'); end
if isempty(Elements),      error('Error - No elements found!'); end

% Is it a single or partitioned domain?
if isempty(PartEntities), single_domain=1; else, single_domain=0; end

% Split data lines into cells (Only in Matlab!)
cells_MF   = splitlines(MeshFormat);
cells_PN   = splitlines(PhysicalNames);
cells_Ent  = splitlines(Entities);
cells_N    = splitlines(Nodes);
cells_E    = splitlines(Elements);
if not(single_domain)
    cells_PEnt = splitlines(PartEntities);
end

%% Define map of costum Boundary types (BC_type)
% BCnames = { 'BCfile','free','wall','outflow',...
%         'imposedPressure','imposedVelocities',...
%         'axisymmetric_y','axisymmetric_x',...
%         'BC_rec','free_rec','wall_rec','outflow_rec',...
%         'imposedPressure_rec','imposedVelocities_rec',...
%         'axisymmetric_y_rec','axisymmetric_x_rec',...
%         'piston_pressure','piston_velocity',...
%         'recordingObject','recObj','piston_stress'};
% BCsolverIds = [0:7,10:20,20,21]; % costume IDs expected in our solver
% BC_type = containers.Map(BCnames,BCsolverIds);

%% Identify critical data within each section:

%********************%
% 1. Read Mesh Format
%********************%
line_data = sscanf(cells_MF{1},'%f %d %d');
info.version   = line_data(1);	% 4.1 is expected
info.file_type = line_data(2);	% 0:ASCII or 1:Binary
info.mode      = line_data(3);	% 1 in binary mode to detect endianness 

% Sanity check
if (info.version ~= 4.1), error('Error - Expected mesh format v4.1'); end
if (info.file_type ~= 0), error('Error - Binary file not allowed'); end

%***********************%
% 2. Read Physical Names
%***********************%
phys = struct('dim',{},'tag',{},'name',{});
numPhysicalNames = sscanf(cells_PN{1},'%d');
for n = 1:numPhysicalNames
   parts = strsplit(cells_PN{n+1});
   phys(n).dim  = str2double(parts{1});
   phys(n).tag  = str2double(parts{2});
   phys(n).name = strrep(parts{3},'"','');
end
mapPhysNames = containers.Map([phys.tag],{phys.name});
info.Dim = max([phys.dim]);

if single_domain
    %***********************%
    % 3. Read Entities
    %***********************%
    l=1; % line counter
    line_data = sscanf(cells_Ent{l},'%d %d %d %d');
    nP = line_data(1);
    nC = line_data(2);
    nS = line_data(3);
    nV = line_data(4);
    l=2; % line counter
    points  = struct('ID',{},'Phys_ID',{});
    curves  = struct('ID',{},'Phys_ID',{});
    surfaces= struct('ID',{},'Phys_ID',{});
    volumes = struct('ID',{},'Phys_ID',{});
    % read points
    if nP>0
        for i = 1:nP
            points(i) = get_entity(cells_Ent{l},'node');
            l = l+1; % update line counter
        end
        point2Phys = containers.Map([points.ID],[points.Phys_ID]);
    end
    % read curves
    if nC>0
        for i = 1:nC
            curves(i) = get_entity(cells_Ent{l},'curve');
            l = l+1; % update line counter
        end
        curve2Phys = containers.Map([curves.ID],[curves.Phys_ID]);
    end
    % read surfaces
    if nS>0
        for i = 1:nS
            surfaces(i) = get_entity(cells_Ent{l},'surface');
            l = l+1; % update line counter
        end
        surf2Phys = containers.Map([surfaces.ID],[surfaces.Phys_ID]);
    end
    % read volumes
    if nV>0
        for i = 1:nV
            volumes(i) = get_entity(cells_Ent{l},'volume');
            l = l+1; % update line counter
        end
        volm2Phys = containers.Map([volumes.ID],[volumes.Phys_ID]);
    end

else
    %******************************%
    % 4. Read Partitioned Entities
    %******************************%
    l=1; info.numPartitions = sscanf(cells_PEnt{l},'%d');
    %l=2; numGhostEntities = sscanf(cells_PEnt{l},'%d'); % not needed
    l=3; line_data = sscanf(cells_PEnt{l},'%d');
    nP = line_data(1);
    nC = line_data(2);
    nS = line_data(3);
    nV = line_data(4);
    l=4; % line counter
    points  = struct('chld_ID',{},'Prnt_ID',{},'Part_ID',{},'Phys_ID',{});
    curves  = struct('chld_ID',{},'Prnt_ID',{},'Part_ID',{},'Phys_ID',{});
    surfaces= struct('chld_ID',{},'Prnt_ID',{},'Part_ID',{},'Phys_ID',{});
    volumes = struct('chld_ID',{},'Prnt_ID',{},'Part_ID',{},'Phys_ID',{});
    % read points
    if nP>0
        for i = 1:nP
            points(i) = get_partitionedEntity(cells_PEnt{l},'node');
            l = l+1; % update line counter
        end
        point2Part = containers.Map([points.chld_ID],{points.Part_ID});
        point2Phys = containers.Map([points.chld_ID],[points.Phys_ID]);
        point2Geom = containers.Map([points.chld_ID],[points.Prnt_ID]);
    end
    % read curves
    if nC>0
        for i = 1:nC
            curves(i) = get_partitionedEntity(cells_PEnt{l},'curve');
            l = l+1; % update line counter
        end
        curve2Part = containers.Map([curves.chld_ID],{curves.Part_ID});
        curve2Phys = containers.Map([curves.chld_ID],[curves.Phys_ID]);
        curve2Geom = containers.Map([curves.chld_ID],[curves.Prnt_ID]);
    end
    % read surfaces
    if nS>0
        for i = 1:nS
            surfaces(i) = get_partitionedEntity(cells_PEnt{l},'surface');
            l = l+1; % update line counter
        end
        surf2Part = containers.Map([surfaces.chld_ID],{surfaces.Part_ID});
        surf2Phys = containers.Map([surfaces.chld_ID],[surfaces.Phys_ID]);
        surf2Geom = containers.Map([surfaces.chld_ID],[surfaces.Prnt_ID]);
    end
    % read volumes
    if nV>0
        for i = 1:nV
            volumes(i) = get_partitionedEntity(cells_PEnt{l},'volume');
            l = l+1; % update line counter
        end
        volm2Part = containers.Map([volumes.chld_ID],{volumes.Part_ID});
        volm2Phys = containers.Map([volumes.chld_ID],[volumes.Phys_ID]);
        volm2Geom = containers.Map([volumes.chld_ID],[volumes.Prnt_ID]);
    end
end

%***********************%
% 5. Read Nodes
%***********************%
l=1; % read first line
line_data = sscanf(cells_N{l},'%d %d %d %d');
numEntityBlocks = line_data(1);
numNodes        = line_data(2);
%minNodeTag      = line_data(3); % not needed
%maxNodeTag      = line_data(4); % not needed
V = get_nodes(cells_N,numEntityBlocks,numNodes);
fprintf('Total vertices found = %g\n',length(V));

%***********************%
% 6. Read Elements
%***********************%
l=1; % read first line
line_data = sscanf(cells_E{l},'%d %d %d %d');
numEntityBlocks = line_data(1);
numElements     = line_data(2);
%minElementsTag  = line_data(3); % not needed
%maxElementsTag  = line_data(4); % not needed

% Allocate space for Elements data
PE.EToV=[]; PE.phys_tag=[]; PE.geom_tag=[]; PE.part_tag=[]; PE.Etype=[];
LE.EToV=[]; LE.phys_tag=[]; LE.geom_tag=[]; LE.part_tag=[]; LE.Etype=[];
SE.EToV=[]; SE.phys_tag=[]; SE.geom_tag=[]; SE.part_tag=[]; SE.Etype=[];
VE.EToV=[]; VE.phys_tag=[]; VE.geom_tag=[]; VE.part_tag=[]; VE.Etype=[];

e0 = 0; % point Element counter
e1 = 0; % Lines Element counter
e2 = 0; % Triangle Element counter
e3 = 0; % Tetrehedron Element counter
for ent = 1:numEntityBlocks
    l = l+1; % update line counter
    line_data = sscanf(cells_E{l},'%d %d %d %d');
    %entityDim = line_data(1); % 0:point, 1:curve, 2:surface, 3:volume
    entityTag = line_data(2); % this is: Entity.ID | Entity.child_ID
    elementType = line_data(3); % 1:line, 2:triangle, 4:tetrahedron, 15:point
    numElementsInBlock = line_data(4);
    %
    for i=1:numElementsInBlock
        l = l+1; % update line counter
        line_data = sscanf(cells_E{l},'%d %d %d %d');
        %elementID = line_data(1); % we use a local numbering instead
        switch elementType % <-- shoudl be entityDim, but we only use 4 element types
            case 1 % Line elements
                e1 = e1 + 1; % update element counter
                LE.Etype(e1,1) = elementType;
                LE.EToV(e1,:) = line_data(2:3);
                LE.phys_tag(e1,1) = curve2Phys(entityTag);
                if not(single_domain)
                    LE.geom_tag(e1,1) = curve2Geom(entityTag);
                    LE.part_tag(e1,1) = curve2Part(entityTag);
                else
                    LE.geom_tag(e1,1) = entityTag;
                end
            case 2 % triangle elements
                e2 = e2 + 1; % update element counter
                SE.Etype(e2,1) = elementType;
                SE.EToV(e2,:) = line_data(2:4);
                SE.phys_tag(e2,1) = surf2Phys(entityTag);
                if not(single_domain)
                    SE.geom_tag(e2,1) = surf2Geom(entityTag);
                    SE.part_tag(e2,1) = surf2Part(entityTag);
                else
                    SE.geom_tag(e2,1) = entityTag;
                end
            case 4 % tetrahedron elements
                e3 = e3 + 1; % update element counter
                VE.Etype(e3,1) = elementType;
                VE.EToV(e3,:) = line_data(2:5);
                VE.phys_tag(e3,1) = volm2Phys(entityTag);
                if not(single_domain)
                    VE.geom_tag(e3,1) = volm2Geom(entityTag);
                    VE.part_tag(e3,1) = volm2Part(entityTag);
                else 
                    VE.geom_tag(e3,1) = entityTag;
                end
            case 15 % Point elements
                e0 = e0 + 1; % update element counter
                PE.Etype(e0,1) = elementType;
                PE.EToV(e0,:) = line_data(2);
                PE.phys_tag(e0,1) = point2Phys(entityTag);
                if not(single_domain)
                    PE.geom_tag(e0,1) = point2Geom(entityTag);
                    PE.part_tag(e0,1) = point2Part(entityTag);
                else
                    PE.geom_tag(e0,1) = entityTag;
                end
            otherwise, error('element not in list');
        end
    end
end
%
fprintf('Total point-elements found = %d\n',e0);
fprintf('Total line-elements found = %d\n',e1);
fprintf('Total surface-elements found = %d\n',e2);
fprintf('Total volume-elements found = %d\n',e3);
% Sanity check
if numElements ~= (e0+e1+e2+e3)
    error('Total number of elements missmatch!'); 
end
%
end % GMSHv4 read function

% Get single entity information:
function entity = get_entity(str_line,type)
    
    vector = str2double(regexp(str_line,'-?[\d.]+(?:e-?\d+)?','match'));

    switch type
        case 'node'
            % 1. get entityTag
            entityTag = vector(1);

            % 3. get entity coordinates % not needed
            %minX(i) = line_data(2); 
            %minY(i) = line_data(3);
            %minZ(i) = line_data(4);

            % 3. get physical tag associated
            numPhysicalTags = vector(5);
            if numPhysicalTags == 0
                physicalTags = NaN;
            else
                physicalTags = zeros(1,numPhysicalTags);
                for j=1:numPhysicalTags
                    physicalTags(j) = vector(5+j);
                end
            end

        otherwise
            % 1. get entityTag
            entityTag = vector(1);
        
            % 2. get entity boxing limits (for visualization) % not needed
            %minX(i) = line_data(2); 
            %minY(i) = line_data(3);
            %minZ(i) = line_data(4);
            %maxX(i) = line_data(5);
            %maxY(i) = line_data(6);
            %maxZ(i) = line_data(7);
            
            % 3. get physical tag associated
            numPhysicalTags = vector(8);
            if numPhysicalTags == 0
                physicalTags = NaN;
            else
                physicalTags = zeros(1,numPhysicalTags);
                for j=1:numPhysicalTags
                    physicalTags(j) = vector(8+j);
                end
            end
        
            % 4. get tags of subentities that define it. % not needed
            %numBoudingEntities = line_data(9+j);
            %entitiesTags = zeros(1,numBoudingEntities);
            %for k=1:numBoudingEntities
            %   entitiesTags(k) = line_data(9+j+k);
            %end
    end
    % output structure:
    entity  = struct('ID',entityTag,'Phys_ID',physicalTags);
end

% Get single partitioned entity information:
function entity = get_partitionedEntity(str_line,type)
    
    vector = str2double(regexp(str_line,'-?[\d.]+(?:e-?\d+)?','match'));

    switch type
        case 'node'
            % 1. get entityTag
            entityTag = vector(1);

            % 2. get parent dimention and tag % no needed
            %parentDim = vector(2);
            parentTag = vector(3);
            
            numPartitionTags = vector(4);
            if numPartitionTags > 1 % <-- mark it as an interface element!
                j=numPartitionTags; partitionTags = Inf;
            else
                j=numPartitionTags; partitionTags = vector(4+j);
            end

            % 3. get entity coordinates % not needed
            %minX = line_data(5+j); 
            %minY = line_data(6+j);
            %minZ = line_data(7+j);

            % 4. get physical tag associated
            numPhysicalTags = vector(8+j);
            if numPhysicalTags == 0 % <-- entity has not physical group!
                physicalTags = NaN;
            else
                physicalTags = zeros(1,numPhysicalTags);
                for k=1:numPhysicalTags
                    physicalTags(k) = vector(8+j+k);
                end
            end

        otherwise
            % 1. get entityTag
            entityTag = vector(1);
        
            % 2. get parent dimention and tag % no needed
            %parentDim = vector(2);
            parentTag = vector(3);
            
            numPartitionTags = vector(4);
            if numPartitionTags > 1 % <-- mark it as an interface element!
                j=numPartitionTags; partitionTags = Inf;
            else
                j=numPartitionTags; partitionTags = vector(4+j);
            end

            % 3. get entity boxing limits (for visualization) % not needed
            %minX(i) = line_data( 5+j); 
            %minY(i) = line_data( 6+j);
            %minZ(i) = line_data( 7+j);
            %maxX(i) = line_data( 8+j); 
            %maxY(i) = line_data( 9+j);
            %maxZ(i) = line_data(10+j);
            
            % 4. get physical tag associated
            numPhysicalTags = vector(11+j);
            if numPhysicalTags == 0 % <-- entity has not physical group!
                physicalTags = NaN;
            else
                physicalTags = zeros(1,numPhysicalTags);
                for k=1:numPhysicalTags
                    physicalTags(k) = vector(11+j+k);
                end
            end
        
            % 5. get tags of subentities that define it. % not needed
            %numBoudingEntities = line_data(12+j+k);
            %entitiesTags = zeros(1,numBoudingEntities);
            %for l=1:numBoudingEntities
            %   entitiesTags(l) = line_data(12+j+k+l);
            %end
    end
    % output structure:
    entity  = struct('chld_ID',entityTag,'Prnt_ID',parentTag,...
                    'Part_ID',partitionTags,'Phys_ID',physicalTags);
end

% Get single partitioned entity information:
function V = get_nodes(cells_N,numNodeBlocks,numNodes)

    % allocate space for nodal data
    V = zeros(numNodes,3); % [x,y,z] 

    l = 1; % this is the parameters line
    % Read nodes blocks:   (can be read in parallel!)
    for ent = 1:numNodeBlocks
        l = l+1; % update line counter
        %
        % Block parameters
        line_data = sscanf(cells_N{l},'%d %d %d %d');
        %entityDim = line_data(1);  % not needed
        %entityTag = line_data(2);  % not needed
        %parametric = line_data(3); % not needed
        numNodesInBlock = line_data(4);
        %
        % Nodes IDs
        nodeTag = zeros(1,numNodesInBlock); % nodeTag
        for i=1:numNodesInBlock
            l = l+1; % update line counter
            nodeTag(i) = sscanf(cells_N{l},'%d');
        end
        %
        % Nodes Coordinates
        for i=1:numNodesInBlock
            l = l+1; % update line counter
            V(nodeTag(i),:) = sscanf(cells_N{l},'%g %g %g'); % [x(i),y(i),z(i)]
        end
    end

end