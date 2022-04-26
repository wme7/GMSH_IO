%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Test GMSH parsers for msh-files written in format version v2.2 and V4.1
%
%      Coded by Manuel A. Diaz @ Pprime | Univ-Poitiers, 2022.01.21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

[V,~,SE,LE,~,~,info] = GMSHparserV2('../meshes/rectangle_v2.msh');
figure(1); subplot(221); viewNodes(V,info);
figure(2); subplot(221); viewLineElements(V,LE.EToV,info);
figure(3); subplot(221); viewSurfaceElements(V,SE.EToV,info);

[V,~,SE,LE,~,~,info] = GMSHparserV4('../meshes/double_piston_v4.msh');
figure(1); subplot(222); viewNodes(V,info);
figure(2); subplot(222); viewLineElements(V,LE.EToV,info);
figure(3); subplot(222); viewSurfaceElements(V,SE.EToV,info);

[V,VE,SE,LE,~,~,info] = GMSHparserV2('../meshes/cuboid_v2.msh');
figure(1); subplot(223); viewNodes(V,info);
figure(2); subplot(223); viewLineElements(V,LE.EToV,info);
figure(3); subplot(223); viewSurfaceElements(V,SE.EToV,info);
figure(4); subplot(221); viewPartVolumes(V,VE.EToV,VE.part_tag,1,info);
figure(4); subplot(222); viewPartVolumes(V,VE.EToV,VE.part_tag,2,info);

[V,VE,SE,LE,PE,~,info] = GMSHparserV4('../meshes/cuboid_v4.msh');
figure(1); subplot(224); viewNodes(V,info);
figure(2); subplot(224); viewLineElements(V,LE.EToV,info);
figure(3); subplot(224); viewSurfaceElements(V,SE.EToV,info);
figure(4); subplot(223); viewPartVolumes(V,VE.EToV,VE.part_tag,1,info);
figure(4); subplot(224); viewPartVolumes(V,VE.EToV,VE.part_tag,2,info);

% Conclusion:
% GMSH format 2.2 is easier to read and to use directly with
% single-partitioned meshes. However, format 4.1 is more suitable for
% handling partitioned domains as it provides the interfacial elements
% required for comunications between partitions. 
%                                                          M.D. 2022.01.21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot/Display mesh elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewNodes(V,info)
    % Plot nodes with global IDs:
    x = V(:,1); y = V(:,2); z = V(:,3);
    h = scatter3(x,y,z,'.r'); 
    hold on
    for i=1:length(V)
        text(x(i),y(i),z(i),num2str(i));
    end
    hold off
    % Print title and axis
    title(sprintf('%d-D GMSH v%g, %d-partitions',info.Dim,info.version,...
        info.numPartitions),'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    if (info.Dim==2),view(2); else 
        zlabel('$z$','Interpreter','latex'); view(3);
    end
    % Use latex font for tick
    set(groot,'defaultAxesTickLabelInterpreter','latex');
end

function viewLineElements(V,EToV,info)
    % Plot surface elements with global IDs:
    x = V(:,1); y = V(:,2); z = V(:,3);
    lx=x(EToV); ly=y(EToV); lz=z(EToV);
    xc=mean(lx,2); yc=mean(ly,2); zc=mean(lz,3);
    hold on
    for i=1:length(EToV)
        plot3(lx(i,:),ly(i,:),lz(i,:),'-r'); 
        text(xc(i),yc(i),zc(i),num2str(i)); 
    end
    hold off
    % Print title and axis
    title(sprintf('%d-D GMSH v%g, %d-partitions',info.Dim,info.version,...
        info.numPartitions),'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    if (info.Dim==2),view(2); else 
        zlabel('$z$','Interpreter','latex'); view(3);
    end
    % Use latex font for tick
    set(groot,'defaultAxesTickLabelInterpreter','latex');
end

function viewSurfaceElements(V,EToV,info)
    % Plot surface elements with global IDs:
    x = V(:,1); y = V(:,2); z = V(:,3);
    trimesh(EToV,x,y,z,'facecolor','none');
    xc=mean(x(EToV),2); yc=mean(y(EToV),2); zc=mean(z(EToV),3);
    hold on
    for i=1:length(EToV)
        %text(xc(i),yc(i),zc(i),num2str(i));
    end
    hold off
    % Print title and axis
    title(sprintf('%d-D GMSH v%g, %d-partitions',info.Dim,info.version,...
        info.numPartitions),'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    if (info.Dim==2),view(2); else 
        zlabel('$z$','Interpreter','latex');
    end
    % Use latex font for tick
    set(groot,'defaultAxesTickLabelInterpreter','latex');
end

function viewPartVolumes(V,EToV,ETags,tag,info)
    % Plot surface elements with global IDs:
    partitioned_EToV=EToV(ETags==tag,:);
    tetramesh(partitioned_EToV,V,'facecolor','w');
    hold off
    % Print title and axis
    title(sprintf('%d-D GMSH v%g, partition %d',info.Dim,info.version,tag),...
        'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    if (info.Dim==2),view(2); else 
        zlabel('$z$','Interpreter','latex');
    end
    % Use latex font for tick
    set(groot,'defaultAxesTickLabelInterpreter','latex');
end