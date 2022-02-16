function drawTetrahedronsFaces(V,EtoV,element,face)

% facemask
nFaces = 4;
f1 = [1,2,3];
f2 = [1,2,4];
f3 = [2,3,4];
f4 = [1,3,4];
facemask = [f1;f2;f3;f4];

% Total number of elements
if element>size(EtoV,1), error('element not found in EtoV'); end
if nFaces~=size(EtoV,2), error('EtoV do not contains valid tetrahedrons'); end

% element vertices
Ex=[V(EtoV(element,1),1), V(EtoV(element,2),1), V(EtoV(element,3),1), V(EtoV(element,4),1)]';
Ey=[V(EtoV(element,1),2), V(EtoV(element,2),2), V(EtoV(element,3),2), V(EtoV(element,4),2)]';
Ez=[V(EtoV(element,1),3), V(EtoV(element,2),3), V(EtoV(element,3),3), V(EtoV(element,4),3)]';

% each element vertices:
v = [Ex,Ey,Ez];

% face data (color)
c = ones(4,3);

scatter3(Ex,Ey,Ez); hold on
%patch('Faces',facemask(face,:),'Vertices',v,'FaceColor','green');
patch('faces',facemask(face,:),'vertices',v,'FaceVertexCData',c,'FaceColor','flat');

hold off; view(3);
axis tight; axis equal;
end