% Calculate spatial contribution to gloab EEG-fMRI correlation and local EEG-fMRI correlation
% SI Table 7 + 8
%
% 2022-08-22 Jonathan Wirsich
[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

for atl=1:1
    %load connectomes
    [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl);
    %init subnetwork masks
    group_size = zeros(1, length(groups));
    for j = 1:length(groups)
        group_size(j) = sum(group_idx==j);
    end

    %load yeo-mapping
    [sameRSN, RSN_labels] = init_yeo(atlases, atl, serialized_path);
    for yeo_it = 1:length(RSN_labels) 
        masks(yeo_it,:) = model.Connectome(sameRSN==yeo_it).vec;
    end
    
    % method to summarize subnetwork 
    % 'spat': average spatial contribution of subnetwork
    % 'corr': FC_EEG-fMRI correlation of subnetwork connections
    methods{1} = 'corr';
    methods{2} = 'spat';
    
    for m = 2:length(methods)
        met = methods{m};
        
        disp(['Compare sumnetworks depedning on metric: ' met])

        maskedtarget_group = zeros(length(groups), length(bands), length(RSN_labels));
        for b = 1:length(bands)
            for j = 1:length(groups) 

                if strcmp(met, 'spat')
                    spatcontrib_group = getSpatialContribution(squeeze(allgroup_mean_EEG(j,b,:))', allgroup_mean_fMRI(j, :));
                end

                for yeo_it = 1:length(RSN_labels) 
                    if strcmp(met, 'corr')
                        maskedtarget_group(j, b , yeo_it) = ...
                            corr(allgroup_mean_fMRI(j,masks(yeo_it,:))', squeeze(allgroup_mean_EEG(j, b,masks(yeo_it,:))));
                        disp([RSN_labels{yeo_it} ' ' int2str(group_size(j)) ' ' groups{j} ...
                            ' -' bands{b} ': ' num2str(corr(allgroup_mean_fMRI(j, masks(yeo_it,:))', squeeze(allgroup_mean_EEG(j, b,masks(yeo_it,:)))))]);
                    elseif strcmp(met, 'spat')
                        maskedtarget_group(j, b, yeo_it) = mean(spatcontrib_group(masks(yeo_it,:)));
                        disp([RSN_labels{yeo_it} ' ' int2str(group_size(j)) ' ' groups{j} ...
                            ' -' bands{b} ': ' num2str(maskedtarget_group(j, b, yeo_it))]);
                    end
                end
            end
        end

        iter = 5000;
        for b = 1:length(bands)
            allsub_fMRI_select = allsub_fMRI(group_idx==1 | group_idx==2,:);
            allsub_EEG_select = squeeze(allsub_EEG(group_idx==1 | group_idx==2,:,:));
            group_size = zeros(2, 1);
            group_size(1) = sum(group_idx==1);
            group_size(2) = sum(group_idx==2);

            %permute
            counts_p = permute_subjectgroups(iter, b, maskedtarget_group([1 2],:,:), met, allsub_fMRI_select, ...
                allsub_EEG_select, group_size, group_idx, masks, 1);
            for mask_it = 1:length(RSN_labels)
               disp(['contr>lTLE (' bands{b} ')- ' RSN_labels{mask_it} ' : ' num2str(counts_p(mask_it)/iter)]) 
            end

            allsub_fMRI_select = allsub_fMRI(group_idx==1 | group_idx==3,:);
            allsub_EEG_select = squeeze(allsub_EEG(group_idx==1 | group_idx==3,:,:));
            group_size = zeros(2, 1);
            group_size(1) = sum(group_idx==1);
            group_size(2) = sum(group_idx==3);
            
            %permute
            counts_p = permute_subjectgroups(iter, b, maskedtarget_group([1 3],:,:), met, allsub_fMRI_select, ...
                allsub_EEG_select, group_size, group_idx, masks, 1);
            for mask_it = 1:length(RSN_labels)
               disp(['contr>rTLE (' bands{b} ')- ' RSN_labels{mask_it} ' : ' num2str(counts_p(mask_it)/iter)]) 
            end

            allsub_fMRI_select = mean(allsub_fMRI(group_idx==3 ,:),1);
            allsub_EEG_select = squeeze(mean(allsub_EEG(group_idx==3,:,:),1));
            group_size = zeros(2, 1);
            group_size(1) = sum(group_idx==3);

        end
    end

end
