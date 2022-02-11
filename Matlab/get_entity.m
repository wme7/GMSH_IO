function [entityTag,physicalTags] = get_entity(str_line,type)
    
    vector = str2double(regexp(str_line,'-?[\d.]+(?:e-?\d+)?','match'));

    switch type
        case 'node'
            % 1. get entityTag
            entityTag = vector(1);

            % 2. get entity coordinates % not needed
            %minX(i) = line_data(2); 
            %minY(i) = line_data(3);
            %minZ(i) = line_data(4);

            % 3. get physical tag associated
            numPhysicalTags = vector(4);
            physicalTags = zeros(1,numPhysicalTags);
            for j=1:numPhysicalTags
                physicalTags(j) = vector(4+j);
            end

        otherwise
            % 1. get entityTag
            entityTag = vector(1);
        
            % 2. get entity limits (for visualization) % not needed
            %minX(i) = line_data(2); 
            %minY(i) = line_data(3);
            %minZ(i) = line_data(4);
            %maxX(i) = line_data(5);
            %maxY(i) = line_data(6);
            %maxZ(i) = line_data(7);
            
            % 3. get physical tag associated
            numPhysicalTags = vector(8);
            physicalTags = zeros(1,numPhysicalTags);
            for j=1:numPhysicalTags
                physicalTags(j) = vector(8+j);
            end
        
            % 4. get tags of subentities that define it. % not needed
            %numBoudingEntities = line_data(9+j);
            %entitiesTags = zeros(1,numBoudingEntities);
            %for k=1:numBoudingEntities
            %   entitiesTags(k) = line_data(9+j+k);
            %end
    end
end