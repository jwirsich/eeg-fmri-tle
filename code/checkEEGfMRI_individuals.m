% check EEG-fMRI for individuals
% SI Table 16 
% 
% 2023-10-13 Jonathan Wirsich

[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

for atl=1:1
    %load connectomes
    [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, ...
        group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl);
    
    hs_idx = zeros(length(group_idx), 1);
    group_meta = [];
    %select metadata
    count = 0;
    for c1 = 1:length(confs)
        for j = 1:length(groups)
            %load meta
            meta = tdfread([serialized_path 'eeg-fmri_' confs{c1} '_' groups{j} '_participants.tsv']);
            meta.participant_id = string(meta.participant_id);
            meta.group = string(meta.group);
            if j == 1
                meta.etiology = strings(length(meta.sex), 1); 
                meta.epilepsy_onset = zeros(length(meta.sex),1); 
                meta.epilepsy_duration = zeros(length(meta.sex),1); 
                meta.spkPerMin = zeros(length(meta.sex),1); 
            else
                meta.etiology = string(meta.etiology);
                meta.spkPerMin = getSpikesPerMinute(confs{c1}, groups{j}, 0);
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
            
            for i = 1:length(meta.age)
                count = count +1;
                %dont look for field in controls (j>1)
                if j > 1 && strcmp(strtrim(meta.etiology(i,:)), 'HS')
                    hs_idx(count) = 1;
                end
            end
        end
    end
        
    for b = 1:5
        for j = [1 2 3]
            vec_size = size(allsub_fMRI, 2);
            groupfMRI = allsub_fMRI(group_idx==j, :);
            groupEEG = squeeze(allsub_EEG(group_idx==j, b, :));

            corrs = zeros(size(groupfMRI, 1), 1);
            for i = 1:size(groupfMRI, 1)
                corrs(i) = corr(groupfMRI(i,:)', groupEEG(i,:)');
            end
            idv_corrEEGfMRI{j} = corrs;
        end
    
        disp(['Mean Ctr: ' num2str(mean(idv_corrEEGfMRI{1})) ' lTLE: ' num2str(mean(idv_corrEEGfMRI{2})) ...
            ' rTLE: ' num2str(mean(idv_corrEEGfMRI{3}))])
        
        [h1, p1] = ttest2(idv_corrEEGfMRI{1}, idv_corrEEGfMRI{2}, 'tail', 'right');
        [h2, p2] = ttest2(idv_corrEEGfMRI{1}, idv_corrEEGfMRI{3}, 'tail', 'left');
        
        %ttest without controlling for covariables
        disp([bands{b} ': ltle<controls: ' num2str(p1)])
        disp([bands{b} ': rtle>controls: ' num2str(p2)])
        
        %control for dataset, sex and age
        corrs_total = corr(allsub_fMRI', squeeze(allsub_EEG(:, b, :))');
        corrs_total = diag(corrs_total);
        
        select_age = group_meta.age(group_idx == 1 | group_idx == 2);
        select_sex = strcmp("F", group_meta.sex(group_idx==1|group_idx==2));
        select_groupconf = conf_idx(group_idx == 1 | group_idx == 2)-1;
        select_corrs = corrs_total(group_idx == 1 | group_idx == 2);
        select_groups = group_idx(group_idx == 1 | group_idx == 2);
        mdl1 = fitlm([select_groups, select_groupconf,select_age, select_sex], select_corrs);
        
        select_age = group_meta.age(group_idx == 1 | group_idx == 3);
        select_sex = strcmp("F", group_meta.sex(group_idx==1|group_idx==3));
        select_groupconf = conf_idx(group_idx == 1 | group_idx == 3)-1;
        select_corrs = corrs_total(group_idx == 1 | group_idx == 3);
        select_groups = group_idx(group_idx == 1 | group_idx == 3);
        mdl2 = fitlm([select_groups, select_groupconf, select_age, select_sex], select_corrs);
        
        disp(['Linear Fit ' bands{b} ': ltle!=controls: ' num2str(mdl1.Coefficients.pValue(2))])
        disp(['Linear Fit ' bands{b} ': rtle!=controls: ' num2str(mdl2.Coefficients.pValue(2))])
    end
    
end