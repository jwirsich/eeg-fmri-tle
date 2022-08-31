classdef Connectome
% CONNECTOME class
% stores a connectome and return matrix or vector (upper triangular) representation
% 22-08-2022 Jonathan Wirsich
    
    properties
        vec
        %this should be probably an atlas object
        regions
    end
    
    methods
        
        % overload constructor: @param vector and @param regions or 
        % @param matrix(regions, regions)
        function obj = Connectome(data, regions)
            if nargin==2
                obj.vec = data;
                %could be 
                obj.regions = regions;
            elseif nargin==1
                obj.vec = obj.getLoopedVec(data);
                obj.regions = size(data, 1);
            end
        end
        
        function mrtx = getMatrix(obj)
            mrtx = zeros(obj.regions);
            count = 1;
            for r1 = 1:obj.regions-1
                for r2 = r1+1:obj.regions
                    mrtx(r1, r2) = obj.vec(count);
                    mrtx(r2, r1) = mrtx(r1, r2);
                    count = count + 1;
                end
            end
        end
    end
    
    methods(Static)
        function vec = getTriuVec(mrtx)
            dim = size(mrtx);
            %convert to triu mask
            mask_ut = triu(true(dim(1),dim(1)),1);
            vec = mrtx(mask_ut);
        end
        
        function vec = getLoopedVec(mrtx)
            %preallocate conn
            dim = size(mrtx);
            
            count = 1;
            %iterate only upper triangle
            for r1 = 1:dim(1)-1
                for r2 = r1+1:dim(1)
                    vec(count) = mrtx(r1,r2);
                    count = count + 1;
                end
            end
        end
    end
    
end

