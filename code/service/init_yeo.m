% Map Desikan atlas to Yeo-7 Atlas labels (exclude subcortical)
% Yeo et al. 2011, JNP
% 
% 22-08-2022 Jonathan Wirsich, Enrico Amico
function [sameRSN, RSN_labels] = init_yeo(atlases, atl, lib_path)

    %Init yeo atlas
    if strcmp(atlases(atl).name, 'desikan')
        tmp = load([lib_path 'aparc_aseg_yeoRS7_68reg_eeg_nosubc_cmfg2dan.mat']);
        regions = 68;
        yeoROIs_eeg = tmp.yeoROIs_eeg;
        yeoOrder_eeg = tmp.yeoOrder_eeg;
    else
        warning('Atlas not implemented')
    end

    RSN_labels = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN'};

    sameRSN = zeros(regions,regions);
    for i=1:length(yeoROIs_eeg)
        for i2=1:length(yeoROIs_eeg)
            if(yeoROIs_eeg(i)==yeoROIs_eeg(i2))
                sameRSN(i,i2) = yeoROIs_eeg(i);
                sameRSN(i2,i) = yeoROIs_eeg(i);
            end
        end
    end
    
end
    