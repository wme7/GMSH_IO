function [V,E,SE,LE,PE,phys_names] = GMSHv2(filename)
% Extract the volume, faces and edges elements contained in a single
% GMSH file in format v2
%
% Coded by Manuel A. Diaz @ d'Alembert.UPMC, 2020.02.15
%
% Example call: GMSHv2('filename.msh')
%
% Output:
%   V: the vertices (nodes coordinates) -- simple array
%   E: volumetric elementes (tetrahedrons) -- structure
%  SE: surface elements (triangles) -- structure
%  LE: linear elements (edges) -- structure
%  PE: point elements (singular vertices) -- structure
%
fID = fopen(filename,'r');
%
phys_names = [];
while (true)
    tline = fgetl(fID);
    if ~ischar(tline); fclose(fID); break; end
    
    if strcmp(tline,'$MeshFormat')
        GMSHformat = parse_format(fID);
    elseif strcmp(tline,'$PhysicalNames')
        phys_names = parse_names(fID,GMSHformat);
    elseif strcmp(tline,'$Nodes')
        V = parse_nodes(fID,GMSHformat);
    elseif strcmp(tline,'$Elements') 
        [E,SE,LE,PE] = parse_elements(fID,GMSHformat);
    end
end
end

function gmshformat = parse_format(fid)
tline = fgetl(fid);
rawformat = sscanf(tline,'%f');
tline = fgetl(fid); % should be EndMeshFormat
gmshformat = rawformat(1);
end

function names = parse_names(fid,gmshformat)
% Line Format:
% physical-dimension physical-number "physical-name"
tline = fgetl(fid);
n_rows = parse_rows(tline,gmshformat);
names = struct('tag',{},'dim',{},'name',{});
for i = 1:n_rows
    tline = fgetl(fid);
    if exist('OCTAVE_VERSION')
        parts = strsplit(tline,' ');
    else
        parts = regexp(tline,' ','split');
    end
    nsz = size(names,2)+1;
    names(nsz).dim = str2double( parts(1) );
    names(nsz).tag = str2double( parts(2) );
    tname = parts(3);
    names(nsz).name = strrep(tname{1},'"','');
end
end

function mat = parse_nodes(fid,gmshformat)
tline = fgetl(fid);
n_rows = parse_rows(tline,gmshformat);
switch floor(gmshformat)
    % Line Format:
    % node-number x-coord y-coord z-coord
    case 2
        mat = fscanf(fid,'%f',[4,n_rows])';
        mat = mat(:,2:4); % let the index be the node id
    case 4
        mat = zeros(n_rows,4);
        while (true)
            tline = fgetl(fid);
            n = sscanf(tline, '%d')';
            n_block = n(4);
            for b = 1:n_block
                tline = fgetl(fid);
                el= sscanf(tline, '%f')';
                mat(el(1),:) = el(1:end);
            end
            if (el(1) == n_rows); break; end % got them all
        end
        tline = fgetl(fid); % get the EndElements
    otherwise; error('cant parse gmsh file of this format');
end
fprintf('Total vertices found = %g\n',n_rows);
end

function n_rows = parse_rows(tline,gmshformat)
n_rows = sscanf(tline,'%d');
switch floor(gmshformat)
    case 2; n_rows = n_rows(1);
    case 4; n_rows = n_rows(2);
    otherwise; error('cant parse gmsh file of this format');
end
end

function [E,SE,LE,PE] = parse_elements(fid,gmshformat)
tline = fgetl(fid);
n_rows = parse_rows(tline,gmshformat);
switch floor(gmshformat)
    case 2; [E,SE,LE,PE] = parse_v2_elements(fid,n_rows);
    otherwise, error('Error: expected GMSH format 2.2 !');
end
end

function [E,SE,LE,PE] = parse_v2_elements(fid,n_rows)
% Line Format:
% elm-number elm-type number-of-tags < tag > ... node-number-list
 E.EToV=[];  E.phys_tag=[];  E.geom_tag=[];  E.part_tag=[];  E.Etype=[];
SE.EToV=[]; SE.phys_tag=[]; SE.geom_tag=[]; SE.part_tag=[]; SE.Etype=[];
LE.EToV=[]; LE.phys_tag=[]; LE.geom_tag=[]; LE.part_tag=[]; LE.Etype=[];
PE.EToV=[]; PE.phys_tag=[]; PE.geom_tag=[]; PE.part_tag=[]; PE.Etype=[];

e0 = 0; % point Element counter
e1 = 0; % Lines Element counter
e2 = 0; % Triangle Element counter
e3 = 0; % Tetrehedron Element counter
for i = 1:n_rows
    tline = fgetl(fid);
    n = sscanf(tline, '%d')';
    %elementID = n(1);
    elementType = n(2);
    numberOfTags = n(3);
    switch elementType
        case 1 % Line elements
            e1 = e1 + 1; % update element counter
            LE.Etype(e1) = elementType;
            LE.EToV(e1,:) = n(3+numberOfTags+1:end);
            if numberOfTags > 0 % get tags if they exist
                tags = n(3+(1:numberOfTags)); % get tags
                % tags(1) : physical entity to which the element belongs
                % tags(2) : elementary number to which the element belongs
                % tags(3) : number of partitions to which the element belongs
                % tags(4) : partition id number
                if length(tags) >= 1
                    LE.phys_tag(e1) = tags(1);
                    if length(tags) >= 2
                        LE.geom_tag(e1) = tags(2);
                        if length(tags) >= 4
                            LE.part_tag(e1) = tags(4);
                        end
                    end
                end
            end
        case 2 % triangle elements
            e2 = e2 + 1; % update element counter
            SE.Etype(e2) = elementType;
            SE.EToV(e2,:) = n(3+numberOfTags+1:end);
            if numberOfTags > 0 % get tags if they exist
                tags = n(3+(1:numberOfTags)); % get tags
                % tags(1) : physical entity to which the element belongs
                % tags(2) : elementary number to which the element belongs
                % tags(3) : number of partitions to which the element belongs
                % tags(4) : partition id number
                if length(tags) >= 1
                    SE.phys_tag(e2) = tags(1);
                    if length(tags) >= 2
                        SE.geom_tag(e2) = tags(2);
                        if length(tags) >= 4
                            SE.part_tag(e2) = tags(4);
                        end
                    end
                end
            end
        case 4 % tetrahedron elements
            e3 = e3 + 1; % update element counter
            E.Etype(e3) = elementType;
            E.EToV(e3,:) = n(3+numberOfTags+1:end);
            if numberOfTags > 0 % get tags if they exist
                tags = n(3+(1:numberOfTags)); % get tags
                % tags(1) : physical entity to which the element belongs
                % tags(2) : elementary number to which the element belongs
                % tags(3) : number of partitions to which the element belongs
                % tags(4) : partition id number
                if length(tags) >= 1
                    E.phys_tag(e3) = tags(1);
                    if length(tags) >= 2
                        E.geom_tag(e3) = tags(2);
                        if length(tags) >= 4
                            E.part_tag(e3) = tags(4);
                        end
                    end
                end
            end
        case 15 % point element
            e0 = e0 + 1; % update element counter
            PE.Etype(e0) = elementType;
            PE.EToV(e0,:) = n(3+numberOfTags+1:end);
            if numberOfTags > 0 % if they exist
                tags = n(3+(1:numberOfTags)); % get tags
                % tags(1) : physical entity to which the element belongs
                % tags(2) : elementary number to which the element belongs
                % tags(3) : number of partitions to which the element belongs
                % tags(4) : partition id number
                if length(tags) >= 1
                    PE.phys_tag(e0) = tags(1);
                    if length(tags) >= 2
                        PE.geom_tag(e0) = tags(2);
                        if length(tags) >= 4
                            PE.part_tag(e0) = tags(4);
                        end
                    end
                end
            end
        otherwise, error('element not set yet!');
    end
end
fprintf('Total point-elements found = %g\n',e0);
fprintf('Total edge/line-elements found = %g\n',e1);
fprintf('Total surface-elements found = %g\n',e2);
fprintf('Total volume-elements found = %g\n',e3);
end