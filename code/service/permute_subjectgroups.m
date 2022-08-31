% permute group_labels of subnetworks
%     targetMetric_group: mean metric to compare to
%     allsub_fMRI_select: selected fMRI (total no of subjects along groups)
%     groupsize: size of each group
%     masks(): mask of subnetworks
%
% 2022-08-22 Jonathan Wirsich
function counts_p = permute_subjectgroups(iter, b, targetMetric_group, targetMetric, allsub_fMRI_select, allsub_EEG_select, group_size, group_idx, masks, group1GEgroup2)

    %vector length of connectome
    conn_vec_length = size(allsub_fMRI_select,2);
    %
    counts_p = zeros(1, size(masks, 1));
    
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
        perm_labels = group_idx(perm_idx);

        %group_idx 1,2 = controls,lTLE
        %TODO this does not take into account rtle vs. ltle idx==2 vs. idx==3
        count_controls = 0;
        
        for i = 1:group_size(1)+group_size(2)
            if perm_labels(i) == 1 && count_controls < group_size(2)
                count_controls = count_controls+1;
                group_fMRI{1} = group_fMRI{1} + allsub_fMRI_select(i,:)/group_size(1);
                group_EEG{1} = group_EEG{1} + squeeze(allsub_EEG_select(i,:,:)/group_size(1));
            elseif perm_labels(i) > 1
                group_fMRI{2} = group_fMRI{2} + allsub_fMRI_select(i,:)/group_size(2);
                group_EEG{2} = group_EEG{2} + squeeze(allsub_EEG_select(i,:,:)/group_size(2));
            end
        end

        if strcmp(targetMetric, 'spat')
            eeg = group_EEG{1};
            perm_metric_eegfmri1 = getSpatialContribution(eeg(b,:)', group_fMRI{1}');
            eeg = group_EEG{2};
            perm_metric_eegfmri2 = getSpatialContribution(eeg(b,:)', group_fMRI{2}');
        end
        
        %break it down to subnetworks (masks)
        for mask_it = 1:size(masks, 1)
            %controls > pat
            if group1GEgroup2
                if strcmp(targetMetric, 'corr')
                    eeg = group_EEG{1};
                    perm_metric_eegfmri1 = corr(eeg(b,masks(mask_it,:))', group_fMRI{1}(masks(mask_it,:))');
                    eeg = group_EEG{2};
                    perm_metric_eegfmri2 = corr(eeg(b,masks(mask_it,:))', group_fMRI{2}(masks(mask_it,:))');
                    diff_patcontro = targetMetric_group(1,b,mask_it) - ...
                        targetMetric_group(2,b,mask_it);
                    if perm_metric_eegfmri1 - perm_metric_eegfmri2 > diff_patcontro
                        counts_p(mask_it) = counts_p(mask_it) + 1;
                    end
                elseif strcmp(targetMetric, 'spat')
                    diff_patcontro = mean(squeeze(targetMetric_group(1,b,mask_it))) - ...
                        mean(squeeze(targetMetric_group(2,b,mask_it)));
                    if (mean(perm_metric_eegfmri1(masks(mask_it,:))) - ...
                        mean(perm_metric_eegfmri2(masks(mask_it,:))) > ...
                        mean(diff_patcontro))
                    
                        counts_p(mask_it) = counts_p(mask_it) +1;
                    end
                end
            else
                if strcmp(targetMetric, 'corr')
                    eeg = group_EEG{1};
                    perm_metric_eegfmri1 = corr(eeg(b,masks(mask_it,:))', group_fMRI{1}(masks(mask_it,:))');
                    eeg = group_EEG{2};
                    perm_metric_eegfmri2 = corr(eeg(b,masks(mask_it,:))', group_fMRI{2}(masks(mask_it,:))');
                    diff_patcontro = targetMetric_group(2,b,mask_it) - ...
                        targetMetric_group(1,b,mask_it);
                    if perm_metric_eegfmri2 - perm_metric_eegfmri1 > diff_patcontro
                        counts_p(mask_it) = counts_p(mask_it) +1;
                    end
                elseif strcmp(targetMetric, 'spat')
                    diff_patcontro = targetMetric_group(2,b,mask_it) - ...
                        targetMetric_group(1,b,mask_it);
                    if (mean(perm_metric_eegfmri2(model.Connectome(sameRSN==yeo_iter).vec)) - ...
                        mean(perm_metric_eegfmri1(model.Connectome(sameRSN==yeo_iter).vec)) > ...
                        mean(diff_patcontro))
                        counts_p(mask_it) = counts_p(mask_it) +1;
                    end
                end
            end
        end

    end

end