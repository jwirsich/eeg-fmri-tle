% Calculate EEG-fMRI spatial FC correlation 
% Data Figure 2 and Supplementary Table 4 + 5
%
% 2022-08-22 Jonathan Wirsich
[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

for atl=1:1
    %load connectomes
    [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl);
    
    disp([atlases(atl).name]);
    for b = 1:length(bands)
        all_corr = corr(mean(allsub_fMRI)', squeeze(mean(allsub_EEG(:,b,:))));
        disp(['All subjects - ' bands{b} ': ' num2str(all_corr)]);
    end
    
    for j = 1:length(groups)
        for b = 1:length(bands)
            all_corr = corr(mean(allsub_fMRI(group_idx==j,:))', squeeze(mean(allsub_EEG(group_idx==j,b,:))));
            disp(['All ' groups{j} ' - ' bands{b} ': ' num2str(all_corr)]);
        end
    end
    
    for co = 1:length(confs)
        for j = 1:length(groups)
            for b = 1:length(bands)
                all_corr = corr(mean(allsub_fMRI(group_idx==j&conf_idx==co,:))', squeeze(mean(allsub_EEG(group_idx==j&conf_idx==co,b,:))));
                disp([confs_label{co} ' ' groups{j} ' - ' bands{b} ': ' num2str(all_corr)]);
            end
        end
    end
    
end
