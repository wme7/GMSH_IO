function drawCommunicationPattern(V,mesh1,mesh2)

% separation
epsilon = 0.05;

% get data from mesh
EtoV1 = mesh1.L_EToV;
comm1 = mesh1.L_CommPattern;

EtoV2 = mesh2.L_EToV;
comm2 = mesh2.L_CommPattern;

% elements in the comm-pattern
K1 = length(comm1);
K2 = length(comm2);
if K1~=K2, error('commPatterns have not the same size'); end

% facemask
facemask(1,:)=[0,1,2]; % Bottom face */
facemask(2,:)=[0,1,3]; % Front face */
facemask(3,:)=[1,2,3]; % Right face */ 
facemask(4,:)=[0,2,3]; % Left face */

facemask(5,:)=[1,2,0]; % Bottom face rotated-1 */
facemask(6,:)=[1,3,0]; % Front face rotated-1 */
facemask(7,:)=[2,3,1]; % Right face rotated-1 */
facemask(8,:)=[2,3,0]; % Left face rotated-1 */

facemask(9 ,:)=[2,0,1]; % Bottom face rotated-2 */
facemask(10,:)=[3,0,1]; % Front face rotated-2 */
facemask(11,:)=[3,1,2]; % Right face rotated-2 */
facemask(12,:)=[3,0,2]; % Left face rotated-2 */

facemask(13,:)=[1,0,2]; % Bottom face reverse */
facemask(14,:)=[1,0,3]; % Front face reverse */
facemask(15,:)=[2,1,3]; % Right face reverse */ 
facemask(16,:)=[2,0,3]; % Left face reverse */

facemask(17,:)=[0,2,1]; % Bottom face reverse-rotated-1 */
facemask(18,:)=[0,3,1]; % Front face reverse-rotated-1 */
facemask(19,:)=[1,3,2]; % Right face reverse-rotated-1 */ 
facemask(20,:)=[0,3,2]; % Left face reverse-rotated-1 */

facemask(21,:)=[2,1,0]; % Bottom face reverse-rotated-2 */
facemask(22,:)=[3,1,0]; % Front face reverse-rotated-2 */
facemask(23,:)=[3,2,1]; % Right face reverse-rotated-2 */ 
facemask(24,:)=[3,2,0]; % Left face reverse-rotated-2 */

% Index correction
facemask=facemask+1;

% face data (color)
c = ones(4,3);

for k=1:K1
    e1 = comm1(k,1);
    f1 = comm1(k,2);
    
    % element vertices
    Ex=[V(EtoV1(e1,1),1), V(EtoV1(e1,2),1), V(EtoV1(e1,3),1), V(EtoV1(e1,4),1)]'+epsilon;
    Ey=[V(EtoV1(e1,1),2), V(EtoV1(e1,2),2), V(EtoV1(e1,3),2), V(EtoV1(e1,4),2)]';
    Ez=[V(EtoV1(e1,1),3), V(EtoV1(e1,2),3), V(EtoV1(e1,3),3), V(EtoV1(e1,4),3)]';

    % each element vertices:
    v1 = [Ex,Ey,Ez];
        
    % Plot the faces
    patch('faces',facemask(f1,:),'vertices',v1,'FaceVertexCData',c,'FaceColor','flat');
    
    e2 = comm1(k,4);
    f2 = comm1(k,5);
    
    % destination element 
    Ex=[V(EtoV2(e2,1),1), V(EtoV2(e2,2),1), V(EtoV2(e2,3),1), V(EtoV2(e2,4),1)]'-epsilon;
    Ey=[V(EtoV2(e2,1),2), V(EtoV2(e2,2),2), V(EtoV2(e2,3),2), V(EtoV2(e2,4),2)]';
    Ez=[V(EtoV2(e2,1),3), V(EtoV2(e2,2),3), V(EtoV2(e2,3),3), V(EtoV2(e2,4),3)]';
    
    % each element vertices:
    v2 = [Ex,Ey,Ez];
    
    % Plot the faces
    patch('faces',facemask(f2,:),'vertices',v2,'FaceVertexCData',c,'FaceColor','flat');
    
    % Plot lines to the nodes
    node1 = [V(EtoV1(e1,facemask(f1,:)),1)+epsilon,V(EtoV1(e1,facemask(f1,:)),2),V(EtoV1(e1,facemask(f1,:)),3)];
    node2 = [V(EtoV2(e2,facemask(f2,:)),1)-epsilon,V(EtoV2(e2,facemask(f2,:)),2),V(EtoV2(e2,facemask(f2,:)),3)];
    
    v = [node1(1,:);node2(1,:)]; pl1=line(v(:,1),v(:,2),v(:,3),'Color','r');
    v = [node1(2,:);node2(2,:)]; pl2=line(v(:,1),v(:,2),v(:,3),'Color','b');
    v = [node1(3,:);node2(3,:)]; pl3=line(v(:,1),v(:,2),v(:,3),'Color','g');
    
    if(k==1), view(3); axis tight; axis equal; end
    
    %waitforbuttonpress;
    pause(1)
    
    pl1.Visible = 'off';
    pl2.Visible = 'off';
    pl3.Visible = 'off';
end

% clf('reset');
% 
% for k=1:K2
%     e1 = comm2(k,1);
%     f1 = comm2(k,2);
%     
%     % element vertices
%     Ex=[V(EtoV2(e1,1),1), V(EtoV2(e1,2),1), V(EtoV2(e1,3),1), V(EtoV2(e1,4),1)]'-epsilon;
%     Ey=[V(EtoV2(e1,1),2), V(EtoV2(e1,2),2), V(EtoV2(e1,3),2), V(EtoV2(e1,4),2)]';
%     Ez=[V(EtoV2(e1,1),3), V(EtoV2(e1,2),3), V(EtoV2(e1,3),3), V(EtoV2(e1,4),3)]';
% 
%     % each element vertices:
%     v1 = [Ex,Ey,Ez];
%         
%     % Plot the faces
%     patch('faces',facemask(f1,:),'vertices',v1,'FaceVertexCData',c,'FaceColor','flat');
%     
%     e2 = comm2(k,4);
%     f2 = comm2(k,5);
%     
%     % destination element 
%     Ex=[V(EtoV1(e2,1),1), V(EtoV1(e2,2),1), V(EtoV1(e2,3),1), V(EtoV1(e2,4),1)]'+epsilon;
%     Ey=[V(EtoV1(e2,1),2), V(EtoV1(e2,2),2), V(EtoV1(e2,3),2), V(EtoV1(e2,4),2)]';
%     Ez=[V(EtoV1(e2,1),3), V(EtoV1(e2,2),3), V(EtoV1(e2,3),3), V(EtoV1(e2,4),3)]';
%     
%     % each element vertices:
%     v2 = [Ex,Ey,Ez];
%     
%     % Plot the faces
%     patch('faces',facemask(f2,:),'vertices',v2,'FaceVertexCData',c,'FaceColor','flat');
%     
%     % Plot lines to the nodes
%     node1 = [V(EtoV2(e1,facemask(f1,:)),1)-epsilon,V(EtoV2(e1,facemask(f1,:)),2),V(EtoV2(e1,facemask(f1,:)),3)];
%     node2 = [V(EtoV1(e2,facemask(f2,:)),1)+epsilon,V(EtoV1(e2,facemask(f2,:)),2),V(EtoV1(e2,facemask(f2,:)),3)];
%     
%     v = [node1(1,:);node2(1,:)]; pl1=line(v(:,1),v(:,2),v(:,3),'Color','r');
%     v = [node1(2,:);node2(2,:)]; pl2=line(v(:,1),v(:,2),v(:,3),'Color','b');
%     v = [node1(3,:);node2(3,:)]; pl3=line(v(:,1),v(:,2),v(:,3),'Color','g');
%     
%     if(k==1), view(3); axis tight; axis equal; end
%     
%     %waitforbuttonpress;
%     pause(1);
%     
%     pl1.Visible = 'off';
%     pl2.Visible = 'off';
%     pl3.Visible = 'off';
% end
hold off;
view(3); axis tight; axis equal;