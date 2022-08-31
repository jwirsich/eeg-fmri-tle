% Define contribution of each node to correlation between EEG and fMRI
% Formula 13 from Colcough 2016 NeuroImage 
%
% Jonathan Wirsich 22-11-2019
function contrib_vec = getSpatialContribution(meanVec1, meanVec2)

    spatialMean1 = mean(meanVec1);
    spatialMean2 = mean(meanVec2);

    denominator = sqrt(sum((meanVec1-spatialMean1).^2));
    z_1 = (meanVec1-spatialMean1)./denominator;

    denominator = sqrt(sum((meanVec2-spatialMean2).^2));
    z_2 = (meanVec2-spatialMean2)./denominator;

    contrib_vec = z_1.*z_2/sum(z_1.*z_2);

end