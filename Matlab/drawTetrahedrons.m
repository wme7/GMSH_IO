function drawTetrahedrons(V,EtoV,EtoBE)

% Total number of elements
K = length(EtoV);

% facemask
f1 = [1,2,3];
f2 = [1,2,4];
f3 = [2,3,4];
f4 = [1,3,4];
f = [f1;f2;f3;f4];

% element vertices
Ex=[V(EtoV(:,1),1), V(EtoV(:,2),1), V(EtoV(:,3),1), V(EtoV(:,4),1)];
Ey=[V(EtoV(:,1),2), V(EtoV(:,2),2), V(EtoV(:,3),2), V(EtoV(:,4),2)];
Ez=[V(EtoV(:,1),3), V(EtoV(:,2),3), V(EtoV(:,3),3), V(EtoV(:,4),3)];

for e=1:K
    % each element vertices:
    v = [Ex(e,:)',Ey(e,:)',Ez(e,:)'];

    % face data (color)
    c = zeros(4,3);
    if EtoBE(e,1)==-1,c(1,:)=[1,1,1]; elseif EtoBE(e,1)==-2;c(1,:)=[0,1,0]; else,c(1,:)=[1,0,0]; end
    if EtoBE(e,2)==-1,c(2,:)=[1,1,1]; elseif EtoBE(e,2)==-2;c(2,:)=[0,1,0]; else,c(2,:)=[1,0,0]; end
    if EtoBE(e,3)==-1,c(3,:)=[1,1,1]; elseif EtoBE(e,3)==-2;c(3,:)=[0,1,0]; else,c(3,:)=[1,0,0]; end
    if EtoBE(e,4)==-1,c(4,:)=[1,1,1]; elseif EtoBE(e,4)==-2;c(4,:)=[0,1,0]; else,c(4,:)=[1,0,0]; end

    patch('faces',f,'vertices',v,'FaceVertexCData',c,'FaceColor','flat'); hold on;
end
hold off; view(3);
axis tight; axis equal;
end