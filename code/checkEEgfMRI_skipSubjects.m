% skip subjects to form subgroups
% SI Table 4, 10, 11, 14
%
% 2023-10-13 Jonathan Wirsich

%select table: options: 'SITable_4', 'SITable_10', 'SITable_11', 'SITable_14_HS', 'SITable_14_noHS'
select_table = 'SITable_4';
[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

for atl=1:1

    %load connectomes
    [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl);
    
    if strcmp(select_table, 'SITable_11') %zscore
        allsub_fMRI = zscore(allsub_fMRI')';
        for band = 1:5
            allsub_EEG(:, band, :) = zscore(squeeze(allsub_EEG(:, band, :))')';
        end 
    end

    %init subnetwork masks
    masks = true(1, size(allsub_fMRI,2));
    mask_labels = cell(1,1);
    mask_labels{1} = 'No Subnetwork Mask';

    group_size = zeros(1, length(groups));
    for j = 1:length(groups)
        group_size(j) = sum(group_idx==j);
    end
    
    hs_idx = zeros(length(group_idx), 1);
    group_meta = [];
    
    %select dataset configuration
    for c1 = 1:length(confs)
        if c1 == 1
            load([serialized_path 'eeg-fmri_256Ch3T_tle_spikes.mat'])
        else
            load([serialized_path 'eeg-fmri_64Ch3T_tle_spikes.mat'])
        end

        for it_subj = 1:length(subj)
           spikecount(it_subj) = length(subj(it_subj).spike_event);
           spikecount_5min(it_subj) = sum(subj(it_subj).spike_event<=150);
        end
        
        for j = 1:length(groups)
            %load meta
            meta = tdfread([serialized_path 'eeg-fmri_' confs{c1} '_' groups{j} '_participants.tsv']);
            meta.participant_id = string(meta.participant_id);
            meta.group = string(meta.group);
            group_idx_before = 0;
            if j == 1
                meta.etiology = strings(length(meta.sex), 1); 
                meta.epilepsy_onset = zeros(length(meta.sex),1); 
                meta.epilepsy_duration = zeros(length(meta.sex),1); 
                meta.IEDmin = zeros(length(meta.sex),1);
            else
                meta.etiology = string(meta.etiology);
                if c1 == 1
                    meta.IEDmin = spikecount_5min(group_idx_before+1:group_idx_before+length(meta.sex))'/5;
                else
                    meta.IEDmin = spikecount(group_idx_before+1:group_idx_before+length(meta.sex))'/20;
                end
                group_idx_before = length(meta.sex);
            end
            
            if 1
                if isempty(group_meta)
                    group_meta = meta;
                else
                    fields = fieldnames(group_meta);
                    for k = 1:numel(fields)
                        aField     = fields{k}; % EDIT: changed to {}
                        merged.(aField) = [group_meta.(aField); meta.(aField)];
                    end
                    group_meta = merged;
                end
            end
            
        end
    end
    
    count = 0;
    for i = 1:length(group_meta.age)
        count = count +1;
        if strcmp(strtrim(group_meta.etiology(i,:)), 'HS')
            hs_idx(count) = 1;
        end
    end
    iter = 5000;

    % select for HS
    group1_select_idx = group_idx==1;
    %select property of patients (e.g HS)
    if strcmp(select_table,  'SITable_14_HS')
        property_select_idx = hs_idx==1;
    elseif strcmp(select_table,  'SITable_14_noHS')
        property_select_idx = hs_idx==0;
    elseif strcmp(select_table,  'SITable_10') %select for spike rate
        property_select_idx = group_meta.IEDmin <= 2;
    else %select all
        property_select_idx = group_meta.IEDmin <= 150;
    end
    
    for it_pat = 2:3 %select epi group to comapre controls against
        tlegr_idx = it_pat;
        if tlegr_idx == 2
            group1GEgroup2 = 1;
        elseif tlegr_idx == 3
            group1GEgroup2 = 0;
        end
        correl_group_select = [1 tlegr_idx];
        group2_select_idx = group_idx==tlegr_idx & property_select_idx;
        select_idx = group1_select_idx | group2_select_idx;

        allsub_fMRI_select = allsub_fMRI(select_idx,:);
        allsub_EEG_select = allsub_EEG(select_idx,:,:);
        group_size(1) = sum(group1_select_idx);
        group_size(2) = sum(group2_select_idx);

        correl_group = zeros(length(groups), length(bands), 1);
        for j = 1:length(groups)
           disp(groups{j})
           for b = 1:length(bands)
               if j==1
                   tmp = allsub_fMRI(group_idx==j,:);
                   tmp = mean(tmp(1:group_size(1),:))';
                   correl_group(j, b, 1) = corr(tmp, ...
                       mean(squeeze(allsub_EEG(group_idx==j,b,:)))');
               else
                   correl_group(j, b, 1) = corr(mean(allsub_fMRI(group_idx==j & property_select_idx,:))', ...
                       mean(squeeze(allsub_EEG(group_idx==j  & property_select_idx,b,:)))');
               end
               disp(bands{b})
               disp(num2str(correl_group(j, b)))
           end
        end

        for b = 1:length(bands)
            %permute
            counts_p = permute_subjectgroups(iter, b, correl_group(correl_group_select,:,:), ...
                'corr', allsub_fMRI_select, allsub_EEG_select, group_size(1:2), ...
                group_idx(select_idx), masks, group1GEgroup2);
            if length(mask_labels) == 1
                if tlegr_idx == 2
                    disp(['Controls > lTLE - ' mask_labels{1} ': ' ...
                        bands{b} ' - ' num2str(counts_p/iter)])
                else
                    disp(['Controls < rTLE - ' mask_labels{1} ': ' ...
                        bands{b} ' - ' num2str(counts_p/iter)])
                end
            else
                for msk_it = 1:length(mask_labels)
                    if tlegr_idx == 2
                        disp(['Controls > lTLE - ' mask_labels{msk_it} ': ' ...
                            bands{b} ' - ' num2str(counts_p(msk_it)/iter)])
                    else
                        disp(['Controls < rTLE - ' mask_labels{msk_it} ': ' ...
                            bands{b} ' - ' num2str(counts_p(msk_it)/iter)])
                    end
                end
            end
        end
    end
    
end
