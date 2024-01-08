% bootstrap average EEG-fMRI correlation and calculate confidence intervals 
% (for Figure 2 and Figure 3)
%
% 2024-01-08 Jonathan Wirsich

[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

for atl=1:1
    %load connectomes
    [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, ...
        group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl);
    
    doDMN = 1; %doDMN = bootsrap DMN correl only (Figure 3)
    if doDMN == 1
        
        [sameRSN, mask_labels] = init_yeo(atlases, atl, serialized_path);
        for yeo_it = 1:length(mask_labels) 
            masks(yeo_it,:) = model.Connectome(sameRSN==yeo_it).vec;
        end
        %get yeo DMN
        allsub_fMRI = allsub_fMRI(:,masks(7,:));
        allsub_EEG = allsub_EEG(:,:,masks(7,:));
    end
    
    for b = 1:length(bands)
        for j = [1 2 3]
            
            test_data = [allsub_fMRI(group_idx==j, :), squeeze(allsub_EEG(group_idx==j, b, :))];
            bootstat = bootstrp(10000, @mean,test_data);

            vec_size = size(allsub_fMRI, 2);
            btstrpd_fMRI = bootstat(:,1:vec_size);
            btstrpd_EEG = bootstat(:,vec_size+1:2*vec_size);

            %corr calcs faster this way than iterating over diagonal
            btsrp_EEGfMRI = corr(btstrpd_fMRI', btstrpd_EEG');
            btsrp_EEGfMRI = diag(btsrp_EEGfMRI);

            % get 95% confidence interval
            p1= prctile(btsrp_EEGfMRI, 97.5);
            p2= prctile(btsrp_EEGfMRI, 2.5);
            mean_btsr = mean(btsrp_EEGfMRI);
            
            diff_p1 = p1-mean_btsr;
            diff_p2 = p2-mean_btsr;
            
            mean_orig = corr(mean(allsub_fMRI(group_idx==j, :))', mean(squeeze(allsub_EEG(group_idx==j, b, :)))');

            disp([bands{b} ' - group ' num2str(j) '; origmean:' num2str(mean_orig) ...
                ' +' num2str(diff_p1) '/' num2str(diff_p2)])
        end
    end
    
end