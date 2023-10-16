% check Eculidian distance vs monomodal EEG or fMRI connectivity
% SI Table 12
%
% 2023-10-13 Jonathan Wirsich

euclDist = getEuclidianDistanceDesikan(); %import Euclidian dstances for Desikan fsaverage
vec_ED = model.Connectome.getLoopedVec(euclDist);

[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

for atl=1:1
    %load connectomes
    [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, ...
        group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl);
    
    %define fMRI as 6th band
    bands{6} = 'fMRI';
    
    datamatrix = zeros(3,5);
    
    for b = 1:6
        for j = 1:3
            if b == 6 %fMRI vs ED
                datamatrix(j, b) = abs(corr(allgroup_mean_fMRI(j,:)', vec_ED'));
            else %EEG vs ED
                datamatrix(j, b) = abs(corr(squeeze(allgroup_mean_EEG(j,b,:)), vec_ED'));
            end
            disp([bands{b} ':' num2str(j) '-' num2str(datamatrix(j, b))]);
        end
    end

    conn_vec_length = size(allsub_fMRI, 2);
    
    for gidx = 2:3
        allsub_fMRI_select = allsub_fMRI(group_idx==1 | group_idx==gidx,:);
        allsub_EEG_select = allsub_EEG(group_idx==1 | group_idx==gidx,:,:);
        group_size(1) = sum(group_idx==1);
        group_size(2) = sum(group_idx==gidx);
        group_idx_select = group_idx(group_idx==1 | group_idx==gidx);

        for b = 1:6
            count_p = 0;
            iter = 5000;
            for it = 1:iter 
                group_fMRI = cell(2,1);
                group_EEG = cell(2,1);
                group_fMRI{1} = zeros(1, conn_vec_length);
                group_fMRI{2} = zeros(1, conn_vec_length);
                group_EEG{1} = zeros(5, conn_vec_length);
                group_EEG{2} = zeros(5, conn_vec_length);

                %ltle beta vs. controls
                perm_idx = randperm(group_size(1)+group_size(2));
                %TODO oversamples groups were only changes within group occur (should happen
                %only once?): offset by requiring result to be strcitly greater than
                %baseline
                perm_labels = group_idx_select(perm_idx);

                %group_idx 1,2 = controls,lTLE
                for i = 1:group_size(1)+group_size(2)
                    if perm_labels(i) == 1
                        group_fMRI{1} = group_fMRI{1} + allsub_fMRI_select(i,:)/group_size(1);
                        group_EEG{1} = group_EEG{1} + squeeze(allsub_EEG_select(i,:,:)/group_size(1));
                    else
                        group_fMRI{2} = group_fMRI{2} + allsub_fMRI_select(i,:)/group_size(2);
                        group_EEG{2} = group_EEG{2} + squeeze(allsub_EEG_select(i,:,:)/group_size(2));
                    end
                end

                if b == 6 %fMRI vs ED
                    corr1 = corr(group_fMRI{1}', vec_ED');
                    corr2 = corr(group_fMRI{2}', vec_ED');
                else %EEG vs ED
                    corr1 = abs(corr(group_EEG{1}(b,:)', vec_ED'));
                    corr2 = abs(corr(group_EEG{2}(b,:)', vec_ED'));
                end

                if gidx == 2
                    if corr1-corr2 > datamatrix(1,b)-datamatrix(gidx,b)
                       count_p = count_p + 1; 
                    end
                elseif gidx == 3
                    if corr2-corr1>datamatrix(gidx,b)-datamatrix(1,b)
                       count_p = count_p + 1; 
                    end
                end

            end

            if gidx == 2
                disp(['ltle<contr -' bands{b} ': p=' num2str(count_p/iter)])
            elseif gidx == 3
                disp(['rtle>contr -' bands{b} ': p=' num2str(count_p/iter)])
            end

        end
    end
end