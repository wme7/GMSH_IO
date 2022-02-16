function [EToE,EToF,EToBE] = buildConnectivities_v4(EToV,BEToV,f1)

nE = size( EToV,1);
nF = size(BEToV,1);
EToE = zeros(nE,1);
EToF = zeros(nE,1);
EToBE= -ones(nE,1);

% Build EtoE and EtoF connectivity
for k = 1:nE
    %    
    V0 = EToV(k,1);
    V1 = EToV(k,2);
    V2 = EToV(k,3);
    V3 = EToV(k,4);
    %
    switch f1
        case 1 % Bottom
            %
            found = 0; % false
            %
            for i = 1:nE
                %
                if (i~=k)
                    %
                    v0 = EToV(i,1);
                    v1 = EToV(i,2);
                    v2 = EToV(i,3);
                    v3 = EToV(i,4);
                    %
                    if((V1 == v0) && (V0 == v1) && (V2 == v2))       % B reversed
                        EToE(k) = i;
                        EToF(k) = 12;
                        found = 1; break;
                    elseif((V1 == v2) && (V0 == v0) && (V2 == v1))   % B reversed rotated 1
                        EToE(k) = i;
                        EToF(k) = 16;
                        found = 1; break;
                    elseif((V1 == v1) && (V0 == v2) && (V2 == v0))   % B reversed rotated 2
                        EToE(k) = i;
                        EToF(k) = 20;
                        found = 1; break;
                    elseif((V1 == v1) && (V0 == v0) && (V2 == v3))   % F straight
                        EToE(k) = i;
                        EToF(k) = 1;
                        found = 1; break;
                    elseif((V1 == v3) && (V0 == v1) && (V2 == v0))   % F reverse rotated 1
                        EToE(k) = i;
                        EToF(k) = 5;
                        found = 1; break;
                    elseif((V1 == v0) && (V0 == v3) && (V2 == v1))   % F reverse rotated 2
                        EToE(k) = i;
                        EToF(k) = 9;
                        found = 1; break;
                    elseif((V1 == v2) && (V0 == v1) && (V2 == v3))   % R straight
                        EToE(k) = i;
                        EToF(k) = 2;
                        found = 1; break;
                    elseif((V1 == v3) && (V0 == v2) && (V2 == v1))   % R reverse rotated 1
                        EToE(k) = i;
                        EToF(k) = 6;
                        found = 1; break;
                    elseif((V1 == v1) && (V0 == v3) && (V2 == v2))   % R reverse rotated 2
                        EToE(k) = i;
                        EToF(k) = 10;
                        found = 1; break;
                    elseif((V1 == v0) && (V0 == v2) && (V2 == v3))   % L straight
                        EToE(k) = i;
                        EToF(k) = 15;
                        found = 1; break;
                    elseif((V1 == v3) && (V0 == v0) && (V2 == v2))   % L reverse rotated 1
                        EToE(k) = i;
                        EToF(k) = 19;
                        found = 1; break;
                    elseif((V1 == v2) && (V0 == v3) && (V2 == v0))   % L reverse rotated 2
                        EToE(k) = i;
                        EToF(k) = 23;
                        found = 1; break;
                    end
                end
            end
            %
            if(not(found))
                EToE(k) =-2; % means: a tetrahedron's face has no connection
                EToF(k) = f1-1;
                for f=1:nF
                    fV1 = BEToV(f,1);
                    fV2 = BEToV(f,2);
                    fV3 = BEToV(f,3);
                    %faceType = BEToV(f,3);
                    
                    if(((V0 == fV1) && (V1 == fV2) && (V2 == fV3)) || ((V0 == fV3) && (V1 == fV1) && (V2 == fV2)) || ((V0 == fV2) && (V1 == fV3) && (V2 == fV1)) || ...
                       ((V0 == fV2) && (V1 == fV1) && (V2 == fV3)) || ((V0 == fV3) && (V1 == fV2) && (V2 == fV1)) || ((V0 == fV1) && (V1 == fV3) && (V2 == fV2)))
                        
                        EToE(k) = k; % means: a tetrahedron's face has connection to a boundary element
                        EToBE(k) = 1; %faceType; % Boundary condition type
                        break;
                    end
                end
            end

        case 2 % Front
            %
            found = 0; % false
            %
            for i = 1:nE
                %
                if (i~=k)
                    %
                    v0 = EToV(i,1);
                    v1 = EToV(i,2);
                    v2 = EToV(i,3);
                    v3 = EToV(i,4);
                    %
                    if((V0 == v0) && (V1 == v1) && (V3 == v2))       % Bottom straight
                        EToE(k) = i;
                        EToF(k) = 0;
                        found = 1; break;
                    elseif((V0 == v2) && (V1 == v0) && (V3 == v1))   % Bottom reverse 1
                        EToE(k) = i;
                        EToF(k) = 8;
                        found = 1; break;
                    elseif((V0 == v1) && (V1 == v2) && (V3 == v0))   % Bottom reverse 2
                        EToE(k) = i;
                        EToF(k) = 4;
                        found = 1; break;
                    elseif((V0 == v1) && (V1 == v0) && (V3 == v3))   % Front straight
                        EToE(k) = i;
                        EToF(k) = 13;
                        found = 1; break;
                    elseif((V0 == v3) && (V1 == v1) && (V3 == v0))   % Front reverse 1
                        EToE(k) = i;
                        EToF(k) = 21;
                        found = 1; break;
                    elseif((V0 == v0) && (V1 == v3) && (V3 == v1))   % Front reverse 2
                        EToE(k) = i;
                        EToF(k) = 17;
                        found = 1; break;
                    elseif((V0 == v2) && (V1 == v1) && (V3 == v3))   % Right straight
                        EToE(k) = i;
                        EToF(k) = 14;
                        found = 1; break;
                    elseif((V0 == v3) && (V1 == v2) && (V3 == v1))   % Right reverse 1
                        EToE(k) = i;
                        EToF(k) = 22;
                        found = 1; break;
                    elseif((V0 == v1) && (V1 == v3) && (V3 == v2))   % Right reverse 2
                        EToE(k) = i;
                        EToF(k) = 18;
                        found = 1; break;
                    elseif((V0 == v0) && (V1 == v2) && (V3 == v3))   % Left straight
                        EToE(k) = i;
                        EToF(k) = 3;
                        found = 1; break;
                    elseif((V0 == v3) && (V1 == v0) && (V3 == v2))   % Left reverse 1
                        EToE(k) = i;
                        EToF(k) = 11;
                        found = 1; break;
                    elseif((V0 == v2) && (V1 == v3) && (V3 == v0))   % Left reverse 2
                        EToE(k) = i;
                        EToF(k) = 7;
                        found = 1; break;
                    end
                end
            end
            %
            if(not(found))
                EToE(k) =-2; % means: a tetrahedron's face has no connection
                EToF(k) = f1-1;
                for f=1:nF
                    fV1 = BEToV(f,1);
                    fV2 = BEToV(f,2);
                    fV3 = BEToV(f,3);
                    %faceType = BEToV(f,3);
                    
                    if(((V0 == fV1) && (V1 == fV2) && (V3 == fV3)) || ((V0 == fV3) && (V1 == fV1) && (V3 == fV2)) || ((V0 == fV2) && (V1 == fV3) && (V3 == fV1)) || ...
                       ((V0 == fV2) && (V1 == fV1) && (V3 == fV3)) || ((V0 == fV3) && (V1 == fV2) && (V3 == fV1)) || ((V0 == fV1) && (V1 == fV3) && (V3 == fV2)))
                        
                        EToE(k) = k; % means: a tetrahedron's face has connection to a boundary element
                        EToBE(k) = 1; %faceType; % Boundary condition type
                        break;
                    end
                end
            end

        case 3 % Right
            %
            found = 0; % false
            %
            for i = 1:nE
                %
                if (i~=k)
                    %
                    v0 = EToV(i,1);
                    v1 = EToV(i,2);
                    v2 = EToV(i,3);
                    v3 = EToV(i,4);
                    %
                    if((V1 == v0) && (V2 == v1) && (V3 == v2))       % Bottom straight
                        EToE(k) = i;
                        EToF(k) = 0;
                        found = 1; break;
                    elseif((V1 == v2) && (V2 == v0) && (V3 == v1))   % Bottom reverse 1
                        EToE(k) = i;
                        EToF(k) = 8;
                        found = 1; break;
                    elseif((V1 == v1) && (V2 == v2) && (V3 == v0))   % Bottom reverse 2
                        EToE(k) = i;
                        EToF(k) = 4;
                        found = 1; break;
                    elseif((V1 == v1) && (V2 == v0) && (V3 == v3))   % Front straight
                        EToE(k) = i;
                        EToF(k) = 13;
                        found = 1; break;
                    elseif((V1 == v3) && (V2 == v1) && (V3 == v0))   % Front reverse 1
                        EToE(k) = i;
                        EToF(k) = 21;
                        found = 1; break;
                    elseif((V1 == v0) && (V2 == v3) && (V3 == v1))   % Front reverse 2
                        EToE(k) = i;
                        EToF(k) = 17;
                        found = 1; break;
                    elseif((V1 == v2) && (V2 == v1) && (V3 == v3))   % Right straight
                        EToE(k) = i;
                        EToF(k) = 14;
                        found = 1; break;
                    elseif((V1 == v3) && (V2 == v2) && (V3 == v1))   % Right reverse 1
                        EToE(k) = i;
                        EToF(k) = 22;
                        found = 1; break;
                    elseif((V1 == v1) && (V2 == v3) && (V3 == v2))   % Right reverse 2
                        EToE(k) = i;
                        EToF(k) = 18;
                        found = 1; break;
                    elseif((V1 == v0) && (V2 == v2) && (V3 == v3))   % Left straight
                        EToE(k) = i;
                        EToF(k) = 3;
                        found = 1; break;
                    elseif((V1 == v3) && (V2 == v0) && (V3 == v2))   % Left reverse 1
                        EToE(k) = i;
                        EToF(k) = 11;
                        found = 1; break;
                    elseif((V1 == v2) && (V2 == v3) && (V3 == v0))   % Left reverse 2
                        EToE(k) = i;
                        EToF(k) = 7;
                        found = 1; break;
                    end
                end
            end
            %
            if(not(found))
                EToE(k) =-2; % means: a tetrahedron's face has no connection
                EToF(k) = f1-1;
                for f=1:nF
                    fV1 = BEToV(f,1);
                    fV2 = BEToV(f,2);
                    fV3 = BEToV(f,3);
                    %faceType = BEToV(f,3);
                    
                    if(((V1 == fV1) && (V2 == fV2) && (V3 == fV3)) || ((V1 == fV3) && (V2 == fV1) && (V3 == fV2)) || ((V1 == fV2) && (V2 == fV3) && (V3 == fV1)) || ...
                       ((V1 == fV2) && (V2 == fV1) && (V3 == fV3)) || ((V1 == fV3) && (V2 == fV2) && (V3 == fV1)) || ((V1 == fV1) && (V2 == fV3) && (V3 == fV2)))
                        
                        EToE(k) = k; % means: a tetrahedron's face has connection to a boundary element
                        EToBE(k) = 1; %faceType; % Boundary condition type
                        break;
                    end
                end
            end

        case 4 % Left
            %
            found = 0; % false
            %
            for i = 1:nE
                %
                if (i~=k)
                    %
                    v0 = EToV(i,1);
                    v1 = EToV(i,2);
                    v2 = EToV(i,3);
                    v3 = EToV(i,4);
                    %
                    if((V2 == v0) && (V0 == v1) && (V3 == v2))       % Bottom straight
                        EToE(k) = i;
                        EToF(k) = 12;
                        found = 1; break;
                    elseif((V2 == v2) && (V0 == v0) && (V3 == v1))   % Bottom reverse 1
                        EToE(k) = i;
                        EToF(k) = 16;
                        found = 1; break;
                    elseif((V2 == v1) && (V0 == v2) && (V3 == v0))   % Bottom reverse 2
                        EToE(k) = i;
                        EToF(k) = 20;
                        found = 1; break;
                    elseif((V2 == v1) && (V0 == v0) && (V3 == v3))   % Front straight
                        EToE(k) = i;
                        EToF(k) = 1;
                        found = 1; break;
                    elseif((V2 == v3) && (V0 == v1) && (V3 == v0))   % Front reverse 1
                        EToE(k) = i;
                        EToF(k) = 5;
                        found = 1; break;
                    elseif((V2 == v0) && (V0 == v3) && (V3 == v1))   % Front reverse 2
                        EToE(k) = i;
                        EToF(k) = 9;
                        found = 1; break;
                    elseif((V2 == v2) && (V0 == v1) && (V3 == v3))   % Right straight
                        EToE(k) = i;
                        EToF(k) = 2;
                        found = 1; break;
                    elseif((V2 == v3) && (V0 == v2) && (V3 == v1))   % Right reverse 1
                        EToE(k) = i;
                        EToF(k) = 6;
                        found = 1; break;
                    elseif((V2 == v1) && (V0 == v3) && (V3 == v2))   % Right reverse 2
                        EToE(k) = i;
                        EToF(k) = 10;
                        found = 1; break;
                    elseif((V2 == v0) && (V0 == v2) && (V3 == v3))   % Left straight
                        EToE(k) = i;
                        EToF(k) = 15;
                        found = 1; break;
                    elseif((V2 == v3) && (V0 == v0) && (V3 == v2))   % Left reverse 1
                        EToE(k) = i;
                        EToF(k) = 19;
                        found = 1; break;
                    elseif((V2 == v2) && (V0 == v3) && (V3 == v0))   % Left reverse 2
                        EToE(k) = i;
                        EToF(k) = 23;
                        found = 1; break;
                    end
                end
            end
            %
            if(not(found))
                EToE(k) =-2; % means: a tetrahedron's face has no connection
                EToF(k) = f1-1;
                for f=1:nF
                    fV1 = BEToV(f,1);
                    fV2 = BEToV(f,2);
                    fV3 = BEToV(f,3);
                    %faceType = BEToV(f,3);
                    
                    if(((V2 == fV1) && (V0 == fV2) && (V3 == fV3)) || ((V2 == fV3) && (V0 == fV1) && (V3 == fV2)) || ((V2 == fV2) && (V0 == fV3) && (V3 == fV1)) || ...
                       ((V2 == fV2) && (V0 == fV1) && (V3 == fV3)) || ((V2 == fV3) && (V0 == fV2) && (V3 == fV1)) || ((V2 == fV1) && (V0 == fV3) && (V3 == fV2)))
                        
                        EToE(k) = k; % means: a tetrahedron's face has connection to a boundary element
                        EToBE(k) = 1; %faceType; % Boundary condition type
                        break;
                    end
                end
            end
    end
end

% Correct EToF for matlab indexing
EToF = EToF+1;

end