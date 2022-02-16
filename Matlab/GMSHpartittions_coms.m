% test GMSHio
clear; close all; clc;

% Read meshfile
[V,E,BE,~,~,physnames] = GMSHparserV2('../meshes/cuboid_v2.msh');

% Get Volume data
 EtoV=E.EToV;   K=length(EtoV);   geoTag=E.geom_tag';   EtoP=E.part_tag';
BEtoV=BE.EToV; KB=length(BEtoV); geoBTag=BE.geom_tag'; BEtoP=BE.part_tag';

% Find elements and faces connectivities
fprintf('Original explicit function: \n');
tic; [EtoE,EtoF,EtoBE] = buildConnectivities(EtoV,BEtoV); toc

% fprintf('Full looped function: \n');
% tic; [EtoE_2,EtoF_2,EtoBE_2] = buildConnectivities_v2(EtoV,BEtoV); toc
% 
% fprintf('Partitioned looped function: \n');
% tic; [EtoE1,EtoF1,EtoBE1] = buildConnectivities_v3(EtoV,BEtoV,1); toc
% tic; [EtoE2,EtoF2,EtoBE2] = buildConnectivities_v3(EtoV,BEtoV,2); toc
% tic; [EtoE3,EtoF3,EtoBE3] = buildConnectivities_v3(EtoV,BEtoV,3); toc
% tic; [EtoE4,EtoF4,EtoBE4] = buildConnectivities_v3(EtoV,BEtoV,4); toc

fprintf('Partitioned explicit function: \n'); 
EtoE_c=zeros(size(EtoE));   EtoF_c=zeros(size(EtoF));   EtoBE_c=zeros(size(EtoBE));
tic; [EtoE_c(:,1),EtoF_c(:,1),EtoBE_c(:,1)] = buildConnectivities_v4(EtoV,BEtoV,1); toc
tic; [EtoE_c(:,2),EtoF_c(:,2),EtoBE_c(:,2)] = buildConnectivities_v4(EtoV,BEtoV,2); toc
tic; [EtoE_c(:,3),EtoF_c(:,3),EtoBE_c(:,3)] = buildConnectivities_v4(EtoV,BEtoV,3); toc
tic; [EtoE_c(:,4),EtoF_c(:,4),EtoBE_c(:,4)] = buildConnectivities_v4(EtoV,BEtoV,4); toc

% Plot the mesh
%figure; drawTetrahedrons(V,EtoV,EtoBE)

%% Find the partitions, their elements and boundaries
 EtoV1=EtoV(EtoP==1,:);    K1=length(EtoV1);
 EtoV2=EtoV(EtoP==2,:);    K2=length(EtoV2);
BEtoV1=BEtoV(BEtoP==1,:); KB1=length(BEtoV1);
BEtoV2=BEtoV(BEtoP==2,:); KB2=length(BEtoV2);

 KP = [ K1,K2 ]; % list of elements/partition
KBP = [KB1,KB2]; % list of boundary elementes/partition

% Find elements and faces connectivities
[EtoE1,EtoF1,EtoBE1] = buildConnectivities(EtoV1,BEtoV1);
[EtoE2,EtoF2,EtoBE2] = buildConnectivities(EtoV2,BEtoV2);

% Plot the partitions
%figure; drawTetrahedrons(V,EtoV1,EtoBE1); title('partition 1');
%figure; drawTetrahedrons(V,EtoV2,EtoBE2); title('partition 2');

% Build mesh-structure for every partition
rank1 = 1;
 % Global elements
mesh1.EToV=EtoV;  
mesh1.EToP=EtoP;  
mesh1.K=K;  
mesh1.KP=KP; 
mesh1.KBP=KBP;
% Local elements
mesh1.L_EToV=EtoV1; 
mesh1.L_BEToV=BEtoV1; 
mesh1.L_EToE=EtoE1; 
mesh1.L_EToF=EtoF1; 
mesh1.L_EToBE=EtoBE1;

rank2 = 2;
% Global elementes
mesh2.EToV=EtoV;  
mesh2.EToP=EtoP;  
mesh2.K=K;  
mesh2.KP=KP; 
mesh2.KBP=KBP; 
% Local elements
mesh2.L_EToV=EtoV2; 
mesh2.L_BEToV=BEtoV2; 
mesh2.L_EToE=EtoE2; 
mesh2.L_EToF=EtoF2; 
mesh2.L_EToBE=EtoBE2;

%% Find elements and faces connectivities between partitions
DEBUG = false;
mesh1 = buildPartitionConnectivities(mesh1,rank1,DEBUG); 
mesh2 = buildPartitionConnectivities(mesh2,rank2,DEBUG); 
%disp(size(mesh1.L_CommPattern));
%disp(size(mesh2.L_CommPattern));

% Plot the partitions and set to green communication elements
%figure; drawTetrahedrons(V,mesh1.L_EToV,mesh1.L_EToBE); title('partition 1');
%figure; drawTetrahedrons(V,mesh2.L_EToV,mesh2.L_EToBE); title('partition 2');

% Sort the communication patterns with respect to their first column
fprintf('commPattern proc1:\n'); disp(sortrows(mesh1.L_CommPattern,[1,3]));
fprintf('commPattern proc2:\n'); disp(sortrows(mesh2.L_CommPattern,[1,3]));

% Draw communication faces
figure; drawCommunicationFaces(V,mesh1); title('partition 1');
figure; drawCommunicationFaces(V,mesh2); title('partition 2');
figure; drawCommunicationPattern(V,mesh1,mesh2);
figure; drawCommunicationPattern(V,mesh2,mesh1);

%%
dir = 'ParadigmS3D_acouBaseTest_2P';
home_path = '~';
binaries_path = '/Depots/paradigms/';
solution_path = ['/Data/DG/',dir,'/'];
INIconfig_path = [pwd,'/tests'];
computing_path = pwd;

%-----------------------------------------
cd([home_path,solution_path]); 
DGmeshInfo_d0.L_KP          = double(h5read('DGmeshInfo_d0.h5','/K'))';
DGmeshInfo_d0.L_KBP         = double(h5read('DGmeshInfo_d0.h5','/KB'))';
DGmeshInfo_d0.L_KBufferi    = double(h5read('DGmeshInfo_d0.h5','/KBuffer'))';
DGmeshInfo_d0.L_EToV        = double(h5read('DGmeshInfo_d0.h5','/EToV'))'+1;
DGmeshInfo_d0.L_EToE        = double(h5read('DGmeshInfo_d0.h5','/EToE'))'+1;
DGmeshInfo_d0.L_EToF        = double(h5read('DGmeshInfo_d0.h5','/EToF'))'+1;
DGmeshInfo_d0.L_EToBE       = double(h5read('DGmeshInfo_d0.h5','/EToBE'))'+1;
DGmeshInfo_d0.L_BEToV       = double(h5read('DGmeshInfo_d0.h5','/BEToV'))'+1;
DGmeshInfo_d0.L_BufferEToV  = double(h5read('DGmeshInfo_d0.h5','/BufferEToV'))'+1;
DGmeshInfo_d0.L_CommPattern = double(h5read('DGmeshInfo_d0.h5','/CommPattern'))'+1;

DGmeshInfo_d1.L_KP          = double(h5read('DGmeshInfo_d1.h5','/K'))';
DGmeshInfo_d1.L_KBP         = double(h5read('DGmeshInfo_d1.h5','/KB'))';
DGmeshInfo_d1.L_KBufferi    = double(h5read('DGmeshInfo_d1.h5','/KBuffer'))';
DGmeshInfo_d1.L_EToV        = double(h5read('DGmeshInfo_d1.h5','/EToV'))'+1;
DGmeshInfo_d1.L_EToE        = double(h5read('DGmeshInfo_d1.h5','/EToE'))'+1;
DGmeshInfo_d1.L_EToF        = double(h5read('DGmeshInfo_d1.h5','/EToF'))'+1;
DGmeshInfo_d1.L_EToBE       = double(h5read('DGmeshInfo_d1.h5','/EToBE'))'+1;
DGmeshInfo_d1.L_BEToV       = double(h5read('DGmeshInfo_d1.h5','/BEToV'))'+1;
DGmeshInfo_d1.L_BufferEToV  = double(h5read('DGmeshInfo_d1.h5','/BufferEToV'))'+1;
DGmeshInfo_d1.L_CommPattern = double(h5read('DGmeshInfo_d1.h5','/CommPattern'))'+1;
cd(computing_path);
%-----------------------------------------

%%
A = (mesh1.KP - DGmeshInfo_d0.L_KP)~=0; nnz(A)
A = (mesh1.KBP - DGmeshInfo_d0.L_KBP)~=0; nnz(A)
A = (mesh1.KBufferi - DGmeshInfo_d0.L_KBufferi)~=0; nnz(A)
A = (mesh1.L_EToV - DGmeshInfo_d0.L_EToV)~=0; nnz(A)
A = (mesh1.L_EToE - DGmeshInfo_d0.L_EToE)~=0; nnz(A)
A = (mesh1.L_EToF - DGmeshInfo_d0.L_EToF)~=0; nnz(A)
A = (mesh1.L_EToBE - DGmeshInfo_d0.L_EToBE)~=0; nnz(A)
%A = (mesh1.L_BEToV - DGmeshInfo_d0.L_BEToV)~=0; nnz(A)
A = (mesh1.L_BufferEToV - DGmeshInfo_d0.L_BufferEToV)~=0; nnz(A)
A = (mesh1.L_CommPattern - DGmeshInfo_d0.L_CommPattern)~=0; nnz(A)

%%
A = (mesh2.KP - DGmeshInfo_d1.L_KP)~=0; nnz(A)
A = (mesh2.KBP - DGmeshInfo_d1.L_KBP)~=0; nnz(A)
A = (mesh2.KBufferi - DGmeshInfo_d1.L_KBufferi)~=0; nnz(A)
A = (mesh2.L_EToV - DGmeshInfo_d1.L_EToV)~=0; nnz(A)
A = (mesh2.L_EToE - DGmeshInfo_d1.L_EToE)~=0; nnz(A)
A = (mesh2.L_EToF - DGmeshInfo_d1.L_EToF)~=0; nnz(A)
A = (mesh2.L_EToBE - DGmeshInfo_d1.L_EToBE)~=0; nnz(A)
%A = (mesh2.L_BEToV - DGmeshInfo_d1.L_BEToV)~=0; nnz(A)
A = (mesh2.L_BufferEToV - DGmeshInfo_d1.L_BufferEToV)~=0; nnz(A)
A = (mesh2.L_CommPattern - DGmeshInfo_d1.L_CommPattern)~=0; nnz(A)
