% Bootstrap Group Means of EEG-fMRI correlation
% Supplementary Table 6
%
% 2022-08-22 Jonathan Wirsich
[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

mdl_names{1} = 'lTLEvsrTLE';
mdl_names{2} = 'controlsVSrTLE';
mdl_names{3} = 'controlsVSlTLE';

for atl=1:1
 
    %load connectomes
    [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl);

    group_size = zeros(1, length(groups));
    for j = 1:length(groups)
        group_size(j) = sum(group_idx==j);
    end
    
    hs_idx = zeros(length(group_idx), 1);
    group_meta = [];
    
    %select
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
   
    iter = 1000;
    for mdlname_idx = 1:length(mdl_names)
        mdl_name = mdl_names{mdlname_idx};
        
        if strcmp(mdl_name, 'lTLEvsrTLE')
            model_labels{1} = 'age';
            model_labels{2} = 'sex';
            model_labels{3} = 'isHS';
            model_labels{4} = 'duration';
            model_labels{5} = 'spkPerMin';
            model_labels{6} = 'ltle0_rtle1';
            model_labels{7} = 'geneva0_marseille1';
        elseif strcmp(mdl_name, 'controlsVSrTLE')
            model_labels{1} = 'age';
            model_labels{2} = 'sex';
            model_labels{3} = 'controls0_rtle1';
            model_labels{4} = 'geneva0_marseille1';
        elseif strcmp(mdl_name, 'controlsVSlTLE')
            model_labels{1} = 'age';
            model_labels{2} = 'sex';
            model_labels{3} = 'controls0_ltle1';
            model_labels{4} = 'geneva0_marseille1';
        elseif strcmp(mdl_name, 'controls')
            model_labels{1} = 'age';
            model_labels{2} = 'sex';
            model_labels{3} = 'geneva0_marseille1';
        end
        p_idxInt = zeros(5, length(model_labels));


        models = cell(5, 1); % one model for each frequency band
        for b = 1:5
            if strcmp(mdl_name, 'lTLEvsrTLE')
                test_data = [allsub_fMRI(group_idx==2|group_idx==3,:), ...
                    squeeze(allsub_EEG(group_idx==2|group_idx==3,b,:)), ...
                    group_meta.age(group_idx==2|group_idx==3), ...
                    strcmp("F", group_meta.sex(group_idx==2|group_idx==3)), ...
                    hs_idx(group_idx==2|group_idx==3), ...
                    group_meta.epilepsy_duration(group_idx==2|group_idx==3), ...
                    group_meta.spkPerMin(group_idx==2|group_idx==3), ...
                    group_idx(group_idx==2|group_idx==3)==3, ...
                    conf_idx(group_idx==2|group_idx==3)==1];
            elseif strcmp(mdl_name, 'controlsVSrTLE')
                test_data = [allsub_fMRI(group_idx==1|group_idx==3,:), ...
                    squeeze(allsub_EEG(group_idx==1|group_idx==3,b,:)), ...
                    group_meta.age(group_idx==1|group_idx==3), ...
                    strcmp("F", group_meta.sex(group_idx==1|group_idx==3)), ...
                    group_idx(group_idx==1|group_idx==3)==3, ...
                    conf_idx(group_idx==1|group_idx==3)==1 ];
            elseif strcmp(mdl_name, 'controlsVSlTLE')
                test_data = [allsub_fMRI(group_idx==1|group_idx==2,:), ...
                    squeeze(allsub_EEG(group_idx==1|group_idx==2,b,:)), ...
                    group_meta.age(group_idx==1|group_idx==2), ...
                    strcmp("F", group_meta.sex(group_idx==1|group_idx==2)), ...
                    group_idx(group_idx==1|group_idx==2)==2, ...
                    conf_idx(group_idx==1|group_idx==2)==1 ];
            elseif strcmp(mdl_name, 'controls')
                test_data = [allsub_fMRI(group_idx==1,:), ...
                    squeeze(allsub_EEG(group_idx==1,b,:)), ...
                    group_meta.age(group_idx==1), ...
                    strcmp("F", group_meta.sex(group_idx==1)), ...
                    conf_idx(group_idx==1)==1 ];
            end

            disp(['reference controls: '  bands{b} ' - ' ...
                num2str(corr(mean(allsub_fMRI(group_idx==1,:),1)', squeeze(mean(allsub_EEG(group_idx==1,b,:),1))))])
            disp(['reference ltle: '  bands{b} ' - ' ...
                num2str(corr(mean(allsub_fMRI(group_idx==2,:),1)', squeeze(mean(allsub_EEG(group_idx==2,b,:),1))))])
            disp(['reference rtle: '  bands{b} ' - ' ...
                num2str(corr(mean(allsub_fMRI(group_idx==3,:),1)', squeeze(mean(allsub_EEG(group_idx==3,b,:),1))))])

            bootstat = bootstrp(iter, @mean,test_data);

            vec_size = size(allsub_fMRI, 2);
            btstrpd_fMRI = bootstat(:,1:vec_size);
            btstrpd_EEG = bootstat(:,vec_size+1:2*vec_size);

            %calculates faster than cycling corrs (only take diagonal)
            btsrp_EEGfMRI = corr(btstrpd_fMRI', btstrpd_EEG');
            btsrp_EEGfMRI = diag(btsrp_EEGfMRI);

            if strcmp(mdl_name, 'lTLEvsrTLE')
                mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                    bootstat(:,vec_size*2+3), bootstat(:,vec_size*2+4) ...
                    bootstat(:,vec_size*2+5), bootstat(:,vec_size*2+6), ...
                    bootstat(:,vec_size*2+7)], btsrp_EEGfMRI);
            elseif strcmp(mdl_name, 'controlsVSrTLE')
                mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                    bootstat(:,vec_size*2+3), bootstat(:,vec_size*2+4)], btsrp_EEGfMRI);
            elseif strcmp(mdl_name, 'controlsVSlTLE')
                mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                    bootstat(:,vec_size*2+3), bootstat(:,vec_size*2+4)], btsrp_EEGfMRI);
            elseif strcmp(mdl_name, 'controls')
                mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                    bootstat(:,vec_size*2+3)], btsrp_EEGfMRI);
            end
            models{b} = mdl;

            %permute labels of interest and put it in a linear model
            %then compare t-values of permuted model to orignial lenar model to
            %extract a p-value
            iter_null = 5000;

            for it_int = 1:length(model_labels)

                count_p = 0;
                idxOfInterest = it_int;

                for it_null = 1:iter_null
                    disp([bands{b} ' idx:' num2str(it_int) ' ' num2str(it_null) '/' num2str(iter_null)]);

                    %TODO not needed we only need lenth of both groups
                    if strcmp(mdl_name, 'lTLEvsrTLE')
                        tmp = group_meta.age(group_idx==2|group_idx==3);
                    elseif strcmp(mdl_name, 'controlsVSrTLE')
                        tmp = group_meta.age(group_idx==1|group_idx==3);
                    elseif strcmp(mdl_name, 'controlsVSlTLE')
                        tmp = group_meta.age(group_idx==1|group_idx==2);
                    elseif strcmp(mdl_name, 'controls')
                        tmp = group_meta.age(group_idx==1);
                    end
                    perm_idx = randperm(length(tmp));

                    %shuffle variable of interest
                    test_data_null = test_data;
                    test_data_null(:, 2278*2+idxOfInterest) = test_data_null(perm_idx, 2278*2+idxOfInterest);

                    bootstat = bootstrp(iter, @mean,test_data_null);

                    vec_size = size(allsub_fMRI, 2);
                    btstrpd_fMRI = bootstat(:,1:vec_size);
                    btstrpd_EEG = bootstat(:,vec_size+1:2*vec_size);

                    %calculates faster than cycling corrs (only take diagonal)
                    btsrp_EEGfMRI = corr(btstrpd_fMRI', btstrpd_EEG');
                    btsrp_EEGfMRI = diag(btsrp_EEGfMRI);

                    if strcmp(mdl_name, 'lTLEvsrTLE')
                        null_mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                            bootstat(:,vec_size*2+3), bootstat(:,vec_size*2+4) ...
                            bootstat(:,vec_size*2+5), bootstat(:,vec_size*2+6), ...
                            bootstat(:,vec_size*2+7)], btsrp_EEGfMRI);
                    elseif strcmp(mdl_name, 'controlsVSrTLE')
                        null_mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                            bootstat(:,vec_size*2+3), bootstat(:,vec_size*2+4)], btsrp_EEGfMRI);
                    elseif strcmp(mdl_name, 'controlsVSlTLE')
                        null_mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                            bootstat(:,vec_size*2+3), bootstat(:,vec_size*2+4)], btsrp_EEGfMRI);
                    elseif strcmp(mdl_name, 'controls')
                        null_mdl = fitlm([bootstat(:,vec_size*2+1), bootstat(:,vec_size*2+2), ...
                            bootstat(:,vec_size*2+3)], btsrp_EEGfMRI);
                    end

                    %get the direction from original model
                    if mdl.Coefficients.tStat(idxOfInterest+1) < 0
                        if null_mdl.Coefficients.tStat(idxOfInterest+1) < mdl.Coefficients.tStat(idxOfInterest+1)
                            count_p = count_p+1;
                        end
                    else
                        if null_mdl.Coefficients.tStat(idxOfInterest+1) > mdl.Coefficients.tStat(idxOfInterest+1)
                            count_p = count_p+1;
                        end
                    end

                end

                disp([bands{b} ' idx:' num2str(it_int) ' p=' num2str(count_p/iter_null)]);
                p_idxInt(b, it_int) = count_p/iter_null;

            end
        end
    
        stat_bootsrapped_model.iter_bootstrap = iter;
        stat_bootsrapped_model.models = models;
        stat_bootsrapped_model.iter_null = iter_null;
        stat_bootsrapped_model.p_idxInt = p_idxInt;
        stat_bootsrapped_model.labels = model_labels;
    
        %serialize
    %     out_path = ['serialized/test_stat_bootsrapped_model_' mdl_name]; 
    %     save(out_path,'stat_bootsrapped_model');
    end

end
