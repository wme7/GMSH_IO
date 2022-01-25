function mesh = GMSHv4(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Extract entities contained in a single GMSH file in format v4. 
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
text = fileread(filename);
%
% Extract strings between:
MeshFormat           = extractBetween(text,'$MeshFormat','$EndMeshFormat');
PhysicalNames        = extractBetween(text,'$PhysicalNames','$EndPhysicalNames');
Entities             = extractBetween(text,'$Entities','$EndEntities');
PartitionedEntities  = extractBetween(text,'$PartitionedEntities','$EndPartitionedEntities');
Nodes                = extractBetween(text,'$Nodes','$EndNodes');
Elements             = extractBetween(text,'$Elements','$EndElements');
Periodic             = extractBetween(text,'$Periodic','$EndPeriodic');
GhostElements        = extractBetween(text,'$GhostElements','$EndGhostElements');
Parametrizations     = extractBetween(text,'$GhostParametrizations','$EndGhostParametrizations');
NodeData             = extractBetween(text,'$NodeData','$EndNodeData');
ElementData          = extractBetween(text,'$ElementData','$EndElementData');
InterpolationMatrices= extractBetween(text,'$InterpolationMatrices','$EndInterpolationMatrices');
%
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
%
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
%
% Identify critical data within each section
%
end % GMSHv4 read function