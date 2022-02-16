function [EToE,EToF,EToBE] = buildConnectivities(EToV,BEToV)

EToE = zeros(size(EToV)); nE=size( EToV,1);
EToF = zeros(size(EToV)); nF=size(BEToV,1);
EToBE= -ones(size(EToV)); 

% Build EtoE and EtoF connectivity
for k = 1:nE
    %
    foundB=0;
    foundF=0;
    foundR=0; 
    foundL=0;
    %    
    V0 = EToV(k,1);
    V1 = EToV(k,2);
    V2 = EToV(k,3);
    V3 = EToV(k,4);
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
            % Bottom
            if((V1 == v0) && (V0 == v1) && (V2 == v2))     % B reversed
                EToE(k,1) = i;
                EToF(k,1) = 12;
                foundB = 1;
            elseif((V1 == v2) && (V0 == v0) && (V2 == v1))   % B reversed rotated 1
                EToE(k,1) = i;
                EToF(k,1) = 16;
                foundB = 1;
            elseif((V1 == v1) && (V0 == v2) && (V2 == v0))   % B reversed rotated 2
                EToE(k,1) = i;
                EToF(k,1) = 20;
                foundB = 1;
            elseif((V1 == v1) && (V0 == v0) && (V2 == v3))   % F straight
                EToE(k,1) = i;
                EToF(k,1) = 1;
                foundB = 1;
            elseif((V1 == v3) && (V0 == v1) && (V2 == v0))   % F reverse rotated 1
                EToE(k,1) = i;
                EToF(k,1) = 5;
                foundB = 1;
            elseif((V1 == v0) && (V0 == v3) && (V2 == v1))   % F reverse rotated 2
                EToE(k,1) = i;
                EToF(k,1) = 9;
                foundB = 1;     
            elseif((V1 == v2) && (V0 == v1) && (V2 == v3))   % R straight
                EToE(k,1) = i;
                EToF(k,1) = 2;
                foundB = 1;
            elseif((V1 == v3) && (V0 == v2) && (V2 == v1))   % R reverse rotated 1
                EToE(k,1) = i;
                EToF(k,1) = 6;
                foundB = 1;
            elseif((V1 == v1) && (V0 == v3) && (V2 == v2))   % R reverse rotated 2
                EToE(k,1) = i;
                EToF(k,1) = 10;
                foundB = 1;
            elseif((V1 == v0) && (V0 == v2) && (V2 == v3))   % L straight
                EToE(k,1) = i;
                EToF(k,1) = 15;
                foundB = 1;
            elseif((V1 == v3) && (V0 == v0) && (V2 == v2))   % L reverse rotated 1
                EToE(k,1) = i;
                EToF(k,1) = 19;
                foundB = 1;
            elseif((V1 == v2) && (V0 == v3) && (V2 == v0))   % L reverse rotated 2
                EToE(k,1) = i;
                EToF(k,1) = 23;
                foundB = 1;
            end

            % Front
            if((V0 == v0) && (V1 == v1) && (V3 == v2))         % Bottom straight
                EToE(k,2) = i;
                EToF(k,2) = 0;
                foundF = 1;
            elseif((V0 == v2) && (V1 == v0) && (V3 == v1))   % Bottom reverse 1
                EToE(k,2) = i;
                EToF(k,2) = 8;
                foundF = 1;                    
            elseif((V0 == v1) && (V1 == v2) && (V3 == v0))   % Bottom reverse 2
                EToE(k,2) = i;
                EToF(k,2) = 4;
                foundF = 1;                                        
            elseif((V0 == v1) && (V1 == v0) && (V3 == v3))   % Front straight
                EToE(k,2) = i;
                EToF(k,2) = 13;
                foundF = 1;
            elseif((V0 == v3) && (V1 == v1) && (V3 == v0))   % Front reverse 1
                EToE(k,2) = i;
                EToF(k,2) = 21;
                foundF = 1;                    
            elseif((V0 == v0) && (V1 == v3) && (V3 == v1))   % Front reverse 2
                EToE(k,2) = i;
                EToF(k,2) = 17;
                foundF = 1;                                        
            elseif((V0 == v2) && (V1 == v1) && (V3 == v3))   % Right straight
                EToE(k,2) = i;
                EToF(k,2) = 14;
                foundF = 1;
            elseif((V0 == v3) && (V1 == v2) && (V3 == v1))   % Right reverse 1
                EToE(k,2) = i;
                EToF(k,2) = 22;
                foundF = 1;                    
            elseif((V0 == v1) && (V1 == v3) && (V3 == v2))   % Right reverse 2
                EToE(k,2) = i;
                EToF(k,2) = 18;
                foundF = 1;                                        
            elseif((V0 == v0) && (V1 == v2) && (V3 == v3))   % Left straight
                EToE(k,2) = i;
                EToF(k,2) = 3;
                foundF = 1;
            elseif((V0 == v3) && (V1 == v0) && (V3 == v2))   % Left reverse 1
                EToE(k,2) = i;
                EToF(k,2) = 11;
                foundF = 1;
            elseif((V0 == v2) && (V1 == v3) && (V3 == v0))   % Left reverse 2
                EToE(k,2) = i;
                EToF(k,2) = 7;
                foundF = 1;
            end

            % Right
            if((V1 == v0) && (V2 == v1) && (V3 == v2))         % Bottom straight
                EToE(k,3) = i;
                EToF(k,3) = 0;
                foundR = 1;
            elseif((V1 == v2) && (V2 == v0) && (V3 == v1))   % Bottom reverse 1
                EToE(k,3) = i;
                EToF(k,3) = 8;
                foundR = 1;
            elseif((V1 == v1) && (V2 == v2) && (V3 == v0))   % Bottom reverse 2
                EToE(k,3) = i;
                EToF(k,3) = 4;
                foundR = 1;
            elseif((V1 == v1) && (V2 == v0) && (V3 == v3))   % Front straight
                EToE(k,3) = i;
                EToF(k,3) = 13;
                foundR = 1;
            elseif((V1 == v3) && (V2 == v1) && (V3 == v0))   % Front reverse 1
                EToE(k,3) = i;
                EToF(k,3) = 21;
                foundR = 1;
            elseif((V1 == v0) && (V2 == v3) && (V3 == v1))   % Front reverse 2
                EToE(k,3) = i;
                EToF(k,3) = 17;
                foundR = 1;
            elseif((V1 == v2) && (V2 == v1) && (V3 == v3))   % Right straight
                EToE(k,3) = i;
                EToF(k,3) = 14;
                foundR = 1;
            elseif((V1 == v3) && (V2 == v2) && (V3 == v1))   % Right reverse 1
                EToE(k,3) = i;
                EToF(k,3) = 22;
                foundR = 1;
            elseif((V1 == v1) && (V2 == v3) && (V3 == v2))   % Right reverse 2
                EToE(k,3) = i;
                EToF(k,3) = 18;
                foundR = 1;
            elseif((V1 == v0) && (V2 == v2) && (V3 == v3))   % Left straight
                EToE(k,3) = i;
                EToF(k,3) = 3;
                foundR = 1;
            elseif((V1 == v3) && (V2 == v0) && (V3 == v2))   % Left reverse 1
                EToE(k,3) = i;
                EToF(k,3) = 11;
                foundR = 1;
            elseif((V1 == v2) && (V2 == v3) && (V3 == v0))   % Left reverse 2
                EToE(k,3) = i;
                EToF(k,3) = 7;
                foundR = 1;
            end

            % Left
            if((V2 == v0) && (V0 == v1) && (V3 == v2))         % Bottom straight
                EToE(k,4) = i;
                EToF(k,4) = 12;
                foundL = 1;
            elseif((V2 == v2) && (V0 == v0) && (V3 == v1))   % Bottom reverse 1
                EToE(k,4) = i;
                EToF(k,4) = 16;
                foundL = 1;
            elseif((V2 == v1) && (V0 == v2) && (V3 == v0))   % Bottom reverse 2
                EToE(k,4) = i;
                EToF(k,4) = 20;
                foundL = 1;
            elseif((V2 == v1) && (V0 == v0) && (V3 == v3))   % Front straight
                EToE(k,4) = i;
                EToF(k,4) = 1;
                foundL = 1;
            elseif((V2 == v3) && (V0 == v1) && (V3 == v0))   % Front reverse 1
                EToE(k,4) = i;
                EToF(k,4) = 5;
                foundL = 1;
            elseif((V2 == v0) && (V0 == v3) && (V3 == v1))   % Front reverse 2
                EToE(k,4) = i;
                EToF(k,4) = 9;
                foundL = 1;
            elseif((V2 == v2) && (V0 == v1) && (V3 == v3))   % Right straight
                EToE(k,4) = i;
                EToF(k,4) = 2;
                foundL = 1;
            elseif((V2 == v3) && (V0 == v2) && (V3 == v1))   % Right reverse 1
                EToE(k,4) = i;
                EToF(k,4) = 6;
                foundL = 1;                        
            elseif((V2 == v1) && (V0 == v3) && (V3 == v2))   % Right reverse 2
                EToE(k,4) = i;
                EToF(k,4) = 10;
                foundL = 1;
            elseif((V2 == v0) && (V0 == v2) && (V3 == v3))   % Left straight
                EToE(k,4) = i;
                EToF(k,4) = 15;
                foundL = 1;
            elseif((V2 == v3) && (V0 == v0) && (V3 == v2))   % Left reverse 1
                EToE(k,4) = i;
                EToF(k,4) = 19;
                foundL = 1;
            elseif((V2 == v2) && (V0 == v3) && (V3 == v0))   % Left reverse 2
                EToE(k,4) = i;
                EToF(k,4) = 23;
                foundL = 1;
            end

            if((foundB == 1) && (foundF == 1) && (foundR == 1) && (foundL == 1))
                break;
            end
        end
    end

    if(foundB == 0)
        EToE(k,1) =-2; % means: a tetrahedron's face has no connection
        EToF(k,1) = 0;
        for f=1:nF
            fV1 = BEToV(f,1);
            fV2 = BEToV(f,2);
            fV3 = BEToV(f,3);
            %faceType = BEToV(f,3);

            if(((V0 == fV1) && (V1 == fV2) && (V2 == fV3)) || ((V0 == fV3) && (V1 == fV1) && (V2 == fV2)) || ((V0 == fV2) && (V1 == fV3) && (V2 == fV1)) || ...
                ((V0 == fV2) && (V1 == fV1) && (V2 == fV3)) || ((V0 == fV3) && (V1 == fV2) && (V2 == fV1)) || ((V0 == fV1) && (V1 == fV3) && (V2 == fV2)))

                EToE(k,1) = k; % means: a tetrahedron's face has connection to a boundary element 
                EToBE(k,1) = 1; %faceType; % Boundary condition type
            end
        end
    end
    if(foundF == 0)
        EToE(k,2) =-2; % means: a tetrahedron's face has no connection
        EToF(k,2) = 1;
        for f=1:nF
            fV1 = BEToV(f,1);
            fV2 = BEToV(f,2);
            fV3 = BEToV(f,3);
            %faceType = BEToV(f,3);

            if(((V0 == fV1) && (V1 == fV2) && (V3 == fV3)) || ((V0 == fV3) && (V1 == fV1) && (V3 == fV2)) || ((V0 == fV2) && (V1 == fV3) && (V3 == fV1)) || ...
                ((V0 == fV2) && (V1 == fV1) && (V3 == fV3)) || ((V0 == fV3) && (V1 == fV2) && (V3 == fV1)) || ((V0 == fV1) && (V1 == fV3) && (V3 == fV2)))

                EToE(k,2) = k; % means: a tetrahedron's face has connection to a boundary element 
                EToBE(k,2) = 1; %faceType; % Boundary condition type
            end
        end
    end
    if(foundR == 0)
        EToE(k,3) =-2; % means: a tetrahedron's face has no connection
        EToF(k,3) = 2;
        for f=1:nF
            fV1 = BEToV(f,1);
            fV2 = BEToV(f,2);
            fV3 = BEToV(f,3);
            %faceType = BEToV(f,3);

            if(((V1 == fV1) && (V2 == fV2) && (V3 == fV3)) || ((V1 == fV3) && (V2 == fV1) && (V3 == fV2)) || ((V1 == fV2) && (V2 == fV3) && (V3 == fV1)) || ...
                ((V1 == fV2) && (V2 == fV1) && (V3 == fV3)) || ((V1 == fV3) && (V2 == fV2) && (V3 == fV1)) || ((V1 == fV1) && (V2 == fV3) && (V3 == fV2)))

                EToE(k,3) = k; % means: a tetrahedron's face has connection to a boundary element 
                EToBE(k,3) = 1; %faceType; % Boundary condition type
            end
        end
    end
    if(foundL == 0)
        EToE(k,4) =-2; % means: a tetrahedron's face has no connection
        EToF(k,4) = 3;
        for f=1:nF
            fV1 = BEToV(f,1);
            fV2 = BEToV(f,2);
            fV3 = BEToV(f,3);
            %faceType = BEToV(f,3);

            if(((V2 == fV1) && (V0 == fV2) && (V3 == fV3)) || ((V2 == fV3) && (V0 == fV1) && (V3 == fV2)) || ((V2 == fV2) && (V0 == fV3) && (V3 == fV1)) || ...
                ((V2 == fV2) && (V0 == fV1) && (V3 == fV3)) || ((V2 == fV3) && (V0 == fV2) && (V3 == fV1)) || ((V2 == fV1) && (V0 == fV3) && (V3 == fV2)))

                EToE(k,4) = k; % means: a tetrahedron's face has connection to a boundary element 
                EToBE(k,4) = 1; %faceType; % Boundary condition type
            end
        end
    end
end

% Correct EToF for matlab indexing
EToF = EToF+1;

end