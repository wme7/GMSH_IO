function drawCommunicationFaces(V,mesh)

% get data from mesh
EtoV = mesh.L_EToV;
comm = mesh.L_CommPattern;

% elements in the comm-pattern
K = length(comm);

% facemask
f1 = [1,2,3];
f2 = [1,2,4];
f3 = [2,3,4];
f4 = [1,3,4];
facemask = [f1;f2;f3;f4];

% face data (color)
c = ones(4,3);

for k=1:K
    element = comm(k,1);
    face    = comm(k,2);
    
    % element vertices
    Ex=[V(EtoV(element,1),1), V(EtoV(element,2),1), V(EtoV(element,3),1), V(EtoV(element,4),1)]';
    Ey=[V(EtoV(element,1),2), V(EtoV(element,2),2), V(EtoV(element,3),2), V(EtoV(element,4),2)]';
    Ez=[V(EtoV(element,1),3), V(EtoV(element,2),3), V(EtoV(element,3),3), V(EtoV(element,4),3)]';

    % each element vertices:
    v = [Ex,Ey,Ez];
    
    % Plot the vertex of the tetrahedron
    scatter3(Ex,Ey,Ez); hold on
    
    % Plot the 
    patch('faces',facemask(face,:),'vertices',v,'FaceVertexCData',c,'FaceColor','flat');
end
hold off;
view(3); axis tight; axis equal;