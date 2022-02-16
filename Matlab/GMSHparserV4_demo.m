% test GMSHio
clear; close all; clc;

% Read meshfile
[V,E,BE,~,~,physnames] = GMSHparserV2('../meshes/cuboid_v2.msh');
%[V,E,BE,~,~,physnames] = GMSHparserV2('../meshes/rectangle_v2.msh');

% Get Volume data
EtoV=E.EToV; nE=size(EtoV,1); geoTag=E.geom_tag'; pTag=E.part_tag';

%---- Visualize 3-d mesh full -----------
%figure; tetramesh(EtoV,V);

%---- Visualize 3-d mesh edges -----------
%figure; trimesh(EtoV,V(:,1),V(:,2),V(:,3),'Edgecolor','b','FaceAlpha','0');
%axis tight; axis equal;

%---- Visualize 3-d mesh partitions -----------
% Ex=[V(EtoV(:,1),1), V(EtoV(:,2),1), V(EtoV(:,3),1)];
% Ey=[V(EtoV(:,1),2), V(EtoV(:,2),2), V(EtoV(:,3),2)];
% Ez=[V(EtoV(:,1),3), V(EtoV(:,2),3), V(EtoV(:,3),3)];
% figure; hold on;
% for e=1:nE
%     if pTag(e)==1
%         color='-c';
%     elseif pTag(e)==2
%         color='-g';
%     else
%         color='-b';
%     end
%     %plot3(Ex(e,[1,2,3,1]),Ey(e,[1,2,3,1]),Ez(e,[1,2,3,1]),color);
%     plot3(Ex(e,[1:end,1]),Ey(e,[1:end,1]),Ez(e,[1:end,1]),color);
% end
% hold off;
% axis tight; axis equal;


% Get Surface data
BEtoV=BE.EToV; nBE=size(BEtoV,1); geoBTag=BE.geom_tag'; pBTag=BE.part_tag';

%---- Visualize 2-d mesh boundaries -----------
Ex=[V(BEtoV(:,1),1), V(BEtoV(:,2),1), V(BEtoV(:,3),1)];
Ey=[V(BEtoV(:,1),2), V(BEtoV(:,2),2), V(BEtoV(:,3),2)];
Ez=[V(BEtoV(:,1),3), V(BEtoV(:,2),3), V(BEtoV(:,3),3)];
figure; hold on;
for e=1:nBE
    if pBTag(e)==1
        color='-c';
    elseif pBTag(e)==2
        color='-g';
    else
        color='-b';
    end
    %plot3(Ex(e,[1,2,3,1]),Ey(e,[1,2,3,1]),Ez(e,[1,2,3,1]),color);
    plot3(Ex(e,[1:end,1]),Ey(e,[1:end,1]),Ez(e,[1:end,1]),color);
end
hold off;
axis tight; axis equal;
