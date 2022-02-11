function [entityTag,partitionTags,physicalTags] = get_partitionedEntity(str_line,type)
    
    vector = str2double(regexp(str_line,'-?[\d.]+(?:e-?\d+)?','match'));

    switch type
        case 'node'
            % 1. get entityTag
            entityTag = vector(1);

            % 2. get parent dimention and tag % no needed
            %parentDim(i) = vector(2);
            %parentTag(i) = vector(3);
            
            numPartitionTags = vector(4);
            partitionTags = zeros(1,numPartitionTags);
            for j=1:numPartitionTags
                partitionTags(j) = vector(4+j);
            end

            % 3. get entity coordinates % not needed
            %minX(i) = line_data(5+j); 
            %minY(i) = line_data(6+j);
            %minZ(i) = line_data(7+j);

            % 4. get physical tag associated
            numPhysicalTags = vector(8+j);
            physicalTags = zeros(1,numPhysicalTags);
            for k=1:numPhysicalTags
                physicalTags(k) = vector(8+j+k);
            end

        otherwise
            % 1. get entityTag
            entityTag = vector(1);
        
            % 2. get parent dimention and tag % no needed
            %parentDim(i) = vector(2);
            %parentTag(i) = vector(3);
            
            numPartitionTags = vector(4);
            partitionTags = zeros(1,numPartitionTags);
            for j=1:numPartitionTags
                partitionTags(j) = vector(4+j);
            end

            % 3. get entity coordinates % not needed
            %minX(i) = line_data( 5+j); 
            %minY(i) = line_data( 6+j);
            %minZ(i) = line_data( 7+j);
            %maxX(i) = line_data( 8+j); 
            %maxY(i) = line_data( 9+j);
            %maxZ(i) = line_data(10+j);
            
            % 4. get physical tag associated
            numPhysicalTags = vector(11+j);
            physicalTags = zeros(1,numPhysicalTags);
            for k=1:numPhysicalTags
                physicalTags(k) = vector(11+j+k);
            end
        
            % 5. get tags of subentities that define it. % not needed
            %numBoudingEntities = line_data(12+j+k);
            %entitiesTags = zeros(1,numBoudingEntities);
            %for l=1:numBoudingEntities
            %   entitiesTags(l) = line_data(12+j+k+l);
            %end
    end
end