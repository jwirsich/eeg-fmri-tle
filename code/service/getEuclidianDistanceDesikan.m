% get ROI center distances from fsaverage Desikan atlas
%
% 2023-10-13 Jonathan Wirsich
function distED = getEuclidianDistanceDesikan()

    regions = 68;
    %load Euclidian distance
    reflect_folder = fileparts(mfilename('fullpath'));
    temp = importdata(fullfile(reflect_folder, '../../data/desi_coord_68.txt'), ' ');

    coord = zeros(length(temp), 3);
    for t = 1:length(temp)
        test = transpose(sscanf(temp{t}, '%f %f %f'));
        coord(t,:) = test;
    end

    distED = zeros(regions);
    for r1 = 1:regions-1
        for r2 = r1+1:regions
            distED(r1, r2) = norm(coord(r1,:)-coord(r2,:));
            distED(r2, r1) = distED(r1, r2);
        end
    end

end