% get subnetwork labels for Desikan atlas
%
% 2023-10-13 Jonathan Wirsich

function [masks, mask_labels] = init_subnetwork_masks()
    %only desikan for now
    regions = 68;
    
    intrahemi_left = zeros(regions);
    for r1 = 1:regions/2
        for r2 = r1+1:regions/2
            intrahemi_left(r1, r2) = 1;
            intrahemi_left(r2, r1) = 1;
        end
    end
    intrahemi_right = zeros(regions);
    for r1 = regions/2+1:regions
        for r2 = r1+1:regions
            intrahemi_right(r1, r2) = 1;
            intrahemi_right(r2, r1) = 1;
        end
    end
    
    masks(1, :) = model.Connectome(intrahemi_left==1).vec;
    mask_labels{1} = 'intrahemisphere_left';
    masks(2, :) = model.Connectome(intrahemi_right==1).vec;
    mask_labels{2} = 'interhemisphere_right';

end