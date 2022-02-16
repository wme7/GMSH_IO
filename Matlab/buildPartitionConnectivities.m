function mesh = buildPartitionConnectivities(mesh,i_rank,DEBUG)

% Partition ID & Total number of Partitions
iP = i_rank; % the rank and the partition number are the same.
nP = length(mesh.KP);

% Elements in the current partition and in the rest of partitions.
iK = mesh.KP(i_rank);
jK = mesh.KP;

% For TET-elements we assume:
Nfaces=4;

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

% Correct indexes for matlab
vn = vn + 1;

% Examing EToE data of the local mesh
[L_E,L_F]=find(mesh.L_EToE==-2);

% Buffer size
K_buffer = length(L_E); % or length(L_F);

% Save to Local mesh parameters
mesh.KBufferi = K_buffer;

% Initialize Buffer elements to vertex array
mesh.L_BufferEToV = zeros(K_buffer,3);

% Find the face vertices for every L_E and accumulate it in BufferToV array
for k=1:K_buffer
    v0=mesh.L_EToV(L_E(k),vn(L_F(k),1));  mesh.L_BufferEToV(k,1)=v0; 
    v1=mesh.L_EToV(L_E(k),vn(L_F(k),2));  mesh.L_BufferEToV(k,2)=v1; 
    v2=mesh.L_EToV(L_E(k),vn(L_F(k),3));  mesh.L_BufferEToV(k,3)=v2;
    if(DEBUG), fprintf('Element|Face : %d | %d, Vertices: ( %d, %d, %d)\n',L_E(k),L_F(k),v0,v1,v2); end
end
%if(DEBUG) std::cout << "BufferEToV" << std::endl; m.L_BufferEToV.print(); end

% Report this connectivities to the local EToBE
for k=1:K_buffer
    mesh.L_EToBE(L_E(k),L_F(k)) = -2; % to denote inter-partition communication
end

% Initialize Buffers of the pattern of communication
BufFedEToE{nP}=[]; % buffer data source
BufFedEToF{nP}=[]; % buffer face source
BufferEToE{nP}=[]; % buffer data destination
BufferEToF{nP}=[]; % buffer face destination
mesh.L_CommPattern = zeros(K_buffer,5); % : (LP_E, LP_F, jP_ID, jP_E, jP_F)

% Build buffer vertex connectivities
for jP=1:nP

    % Clear temporal containers
    Elems=[]; Faces=[]; Efrom=[]; Ffrom=[];

    if (iP~=jP)

        %if(DEBUG) {std::cout << "searching connections in partition " << jP << std::endl;} end

        % Isolate the j-partition of the global EToV
        j_EToV = zeros(jK(jP),Nfaces);
        
        n=1;
        for k=1:mesh.K
            if (mesh.EToP(k)==jP)
                j_EToV(n,1)=mesh.EToV(k,1);
                j_EToV(n,2)=mesh.EToV(k,2);
                j_EToV(n,3)=mesh.EToV(k,3);
                j_EToV(n,4)=mesh.EToV(k,4);
                n=n+1;
            end
        end
        %if(DEBUG) {std::cout << "j_EToV " << std::endl; j_EToV.print();} end

        % Initialize the variables for brute_force matching
        E=jK(jP); 

        % Search for the connecting element and its face
        for k=1:K_buffer

            % load vertices of buffer element
            V0 = mesh.L_BufferEToV(k,1);
            V1 = mesh.L_BufferEToV(k,2);
            V2 = mesh.L_BufferEToV(k,3);

            % brute_force matching faces
            e=1; % element counter

            while (e <= E)                
                % Examing any match with the bottom face
                i = 1;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e)); 
                    Faces = [Faces,i]; %.push_back(0); % the bottom face (facemask id)
                    break;
                end

                % Examing any match with the front face
                i = 2;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(1); % the front face (facemask id)
                    break;
                end

                % Examing any match with the right face
                i = 3;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(2); % the right face (facemask id)
                    break;
                end

                % Examing any match with the left face
                i = 4;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(3); % the left face (facemask id)
                    break;
                end

                % Examing any match with the rotated-1 bottom face
                i = 5;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(4); % the bottom face (facemask id)
                    break;
                end

                % Examing any match with the rotated-1 front face
                i = 6;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(5); % the front face (facemask id)
                    break;
                end

                % Examing any match with the rotated-1 right face
                i = 7;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(6); % the right face (facemask id)
                    break;
                end

                % Examing any match with the rotated-1 left face
                i = 8;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(7); % the left face (facemask id)
                    break;
                end

                % Examing any match with the rotated-2 bottom face
                i = 9;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(8); % the bottom face (facemask id)
                    break;
                end

                % Examing any match with the rotated-2 front face
                i = 10;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(9); % the front face (facemask id)
                    break;
                end

                % Examing any match with the rotated-2 right face
                i = 11;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(10); % the right face (facemask id)
                    break;
                end

                % Examing any match with the rotated-2 left face
                i = 12;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(11); % the left face (facemask id)
                    break;
                end

                % Examing any match with the reversed-bottom face
                i = 13;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(12); % the bottom face (facemask id)
                    break;
                end

                % Examing any match with the reversed-front face
                i = 14;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(13); % the front face (facemask id)
                    break;
                end

                % Examing any match with the reversed-right face
                i = 15;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(14); % the right face (facemask id)
                    break;
                end

                % Examing any match with the reversed-left face
                i = 16;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(15); % the left face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-1 bottom face
                i = 17;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(16); % the bottom face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-1 front face
                i = 18;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e)); 
                    Faces = [Faces,i]; %.push_back(17); % the front face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-1 right face
                i = 19;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(18); % the right face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-1 left face
                i = 20;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(19); % the left face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-2 bottom face
                i = 21;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(20); % the bottom face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-2 front face
                i = 22;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %s.push_back(21); % the front face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-2 right face
                i = 23;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(22); % the right face (facemask id)
                    break;
                end

                % Examing any match with the reversed-rotated-2 left face
                i = 24;
                eV0 = j_EToV(e,vn(i,1));
                eV1 = j_EToV(e,vn(i,2));
                eV2 = j_EToV(e,vn(i,3));

                if((V0 == eV0) && (V1 == eV1) && (V2 == eV2))
                    if(DEBUG), fprintf('Local_partition_%d_E: %d with face %d connects with partition_%d_E: %d through its face %d\n',iP,L_E(k),L_F(k),jP,e,i); end
                    Efrom = [Efrom,L_E(k)]; %.push_back(L_E(k));
                    Ffrom = [Ffrom,L_F(k)]; %.push_back(L_F(k));
                    Elems = [Elems,e]; %.push_back(int(e));
                    Faces = [Faces,i]; %.push_back(23); % the left face (facemask id)
                    break;
                end

                % No match found -> go to next element
                e=e+1;
            end
        end
    end
    BufFedEToE{jP}=Efrom; %.push_back(Efrom);
    BufFedEToF{jP}=Ffrom; %.push_back(Ffrom);
    BufferEToE{jP}=Elems; %.push_back(Elems);
    BufferEToF{jP}=Faces; %.push_back(Faces);
    % Verify 
    if(DEBUG), end%{std::cout << "BufFedEToE" << std::endl; for (auto v:BufFedEToE) {for (int x:v) std::cout << x << ' '; std::cout << std::endl;}}
    if(DEBUG), end%{std::cout << "BufFedEToF" << std::endl; for (auto v:BufFedEToF) {for (int x:v) std::cout << x << ' '; std::cout << std::endl;}}
    if(DEBUG), end%{std::cout << "BufferEToE" << std::endl; for (auto v:BufferEToE) {for (int x:v) std::cout << x << ' '; std::cout << std::endl;}}
    if(DEBUG), end%{std::cout << "BufferEToF" << std::endl; for (auto v:BufferEToF) {for (int x:v) std::cout << x << ' '; std::cout << std::endl;}}
end

% Flattening Info and Building Connectivity patterns
n=1;
for i=1:length(BufferEToE)
    for j=1:length(BufferEToE{i})
        % mesh.L_CommPattern(n,0) = iP; % LP_ID : Local_partition
        mesh.L_CommPattern(n,1) = BufFedEToE{i}(j); % LP_E : Local_partition_Element
        mesh.L_CommPattern(n,2) = BufFedEToF{i}(j); % LP_F : Local_partition_Face
        mesh.L_CommPattern(n,3) = i; % jP_ID : Destination_partition
        mesh.L_CommPattern(n,4) = BufferEToE{i}(j); % jP_E : Destination_partition_Element
        mesh.L_CommPattern(n,5) = BufferEToF{i}(j); % jP_F : Destination_partition_Face
        n=n+1;
    end
end
%if(DEBUG) {std::cout << "CommPattern" << std::endl; m.L_CommPattern.print();} end
mesh.L_CommPattern = sortrows(mesh.L_CommPattern,[1,3]); % sort by origin element and Partition
end