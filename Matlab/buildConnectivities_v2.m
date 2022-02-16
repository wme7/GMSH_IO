function [EToE,EToF,EToBE] = buildConnectivities_v2(EToV,BEToV)

EToE = zeros(size(EToV)); nE=size( EToV,1);
EToF = zeros(size(EToV)); nF=size(BEToV,1);
EToBE= -ones(size(EToV)); 

% For TET-elements we assume:
Nfaces=4; Nfmasks=24; Nbfmasks=6;

% Facemask
vn(1,:)=[0,1,2]; % Bottom face */
vn(2,:)=[0,1,3]; % Front face */
vn(3,:)=[1,2,3]; % Right face */ 
vn(4,:)=[0,2,3]; % Left face */

vn(5,:)=[1,2,0]; % Bottom face rotated-1 */
vn(6,:)=[1,3,0]; % Front face rotated-1 */
vn(7,:)=[2,3,1]; % Right face rotated-1 */
vn(8,:)=[2,3,0]; % Left face rotated-1 */

vn(9 ,:)=[2,0,1]; % Bottom face rotated-2 */
vn(10,:)=[3,0,1]; % Front face rotated-2 */
vn(11,:)=[3,1,2]; % Right face rotated-2 */
vn(12,:)=[3,0,2]; % Left face rotated-2 */

vn(13,:)=[1,0,2]; % Bottom face reverse */
vn(14,:)=[1,0,3]; % Front face reverse */
vn(15,:)=[2,1,3]; % Right face reverse */ 
vn(16,:)=[2,0,3]; % Left face reverse */

vn(17,:)=[0,2,1]; % Bottom face reverse-rotated-1 */
vn(18,:)=[0,3,1]; % Front face reverse-rotated-1 */
vn(19,:)=[1,3,2]; % Right face reverse-rotated-1 */ 
vn(20,:)=[0,3,2]; % Left face reverse-rotated-1 */

vn(21,:)=[2,1,0]; % Bottom face reverse-rotated-2 */
vn(22,:)=[3,1,0]; % Front face reverse-rotated-2 */
vn(23,:)=[3,2,1]; % Right face reverse-rotated-2 */ 
vn(24,:)=[3,2,0]; % Left face reverse-rotated-2 */

vb(1,:)=vn( 1,:); % Bottom face */
vb(2,:)=vn( 5,:); % Bottom face rotated-1 */
vb(3,:)=vn( 9,:); % Bottom face rotated-2 */
vb(4,:)=vn(13,:); % Bottom face reverse */
vb(5,:)=vn(17,:); % Bottom face reverse-rotated-1 */
vb(6,:)=vn(21,:); % Bottom face reverse-rotated-2 */

% Correct indexes for matlab
vn = vn + 1;
vb = vb + 1;

% Build EtoE and EtoF connectivity for k-element
for k = 1:nE
    %
    % Search for each face of k-element
    for f1 = 1:Nfaces
        %
        found = false;
        %
        V0 = EToV(k,vn(f1,1));
        V1 = EToV(k,vn(f1,2));
        V2 = EToV(k,vn(f1,3));
        %
        % Search all element's faces
        for i = 1:nE
            %
            if (i~=k)
                %
                for f2 = 1:Nfmasks
                    %
                    v0 = EToV(i,vn(f2,1));
                    v1 = EToV(i,vn(f2,2));
                    v2 = EToV(i,vn(f2,3));
                    %
                    if((V0 == v0) && (V1 == v1) && (V2 == v2))
                        EToE(k,f1) = i;
                        EToF(k,f1) = f2-1; % correct for matlab's indexing
                        found = true;
                        break; % go to endOfLoop
                    end
                end
            end
            if (found), break; end
        end
        % End of search all element's faces
        %
        % Search on the boundary (face) elements
        if(~found)
            EToE(k,f1) =-2; % means: a tetrahedron's face has no connection
            EToF(k,f1) =f1-1; % correct for matlab's indexing
            for f=1:nF
                %
                for i=1:Nbfmasks
                    %
                    fV1 = BEToV(f,vb(i,1));
                    fV2 = BEToV(f,vb(i,2));
                    fV3 = BEToV(f,vb(i,3));
                    %faceType = BEToV(f,3);
                    if ((V0 == fV1) && (V1 == fV2) && (V2 == fV3))
                        EToE(k,f1) = k; % means: a tetrahedron's face has connection to a boundary element 
                        EToBE(k,f1) = 1; % faceType; % Boundary condition type
                        found = true;
                        break; % go to endOfLoop
                    end
                end
                if (found), break; end
            end
        end
        % EndOfLoop
    end
end

% Correct EToF for matlab indexing
EToF = EToF+1;

end