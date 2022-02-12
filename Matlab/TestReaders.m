%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Test GMSH readers for format version v2.2 (Legacy) and V4.1
%
%      Coded by Manuel A. Diaz @ d'Alembert | UPMC, 2020.02.15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

[V,~,SE,LE,~,~,info] = GMSHv2('../meshes/rectangle_v2.msh');
viewInfo(info); % Display info
figure(1); subplot(221); viewNodes(V,info);
figure(2); subplot(221); viewLineElements(V,LE.EToV,info);
figure(3); subplot(221); viewSurfaceElements(V,SE.EToV,info);

[V,~,SE,LE,~,~,info] = GMSHv4('../meshes/rectangle_v4.msh');
viewInfo(info); % Display info
figure(1); subplot(222); viewNodes(V,info);
figure(2); subplot(222); viewLineElements(V,LE.EToV,info);
figure(3); subplot(222); viewSurfaceElements(V,SE.EToV,info);

[V,VE,SE,LE,~,~,info] = GMSHv2('../meshes/cuboid_v2.msh');
viewInfo(info); % Display info
figure(1); subplot(223); viewNodes(V,info);
figure(2); subplot(223); viewLineElements(V,LE.EToV,info);
figure(3); subplot(223); viewSurfaceElements(V,SE.EToV,info);
figure(4); subplot(221); viewPartVolumes(V,VE.EToV,VE.part_tag,1,info);
figure(4); subplot(222); viewPartVolumes(V,VE.EToV,VE.part_tag,2,info);

[V,VE,SE,LE,PE,~,info] = GMSHv4('../meshes/cuboid_v4.msh');
viewInfo(info); % Display info
figure(1); subplot(224); viewNodes(V,info);
figure(2); subplot(224); viewLineElements(V,LE.EToV,info);
figure(3); subplot(224); viewSurfaceElements(V,SE.EToV,info);
figure(4); subplot(223); viewPartVolumes(V,VE.EToV,VE.part_tag,1,info);
figure(4); subplot(224); viewPartVolumes(V,VE.EToV,VE.part_tag,2,info);

% Conclusion:
% GMSH format 2.2 is easier to read and use for single-partitioned meshes.
% However, format 4.1 is more suitable for handling partitioned domains.
%                                                          M.D. 2022.02.08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot/Display mesh data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewInfo(info)
    % Display mesh format information:
    fprintf('Mesh version %g, Binary %d, endian %d\n',...
            info.version,info.file_type,info.mode);
end

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