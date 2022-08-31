% load all EEG and fMRI datasets
%
% 2022-08-22 Jonathan Wirsich
function [allgroup_mean_fMRI, allgroup_mean_EEG, allsub_fMRI, allsub_EEG, groups, group_idx, conf_idx, bands, atlases] = loadDatasets4eeg_fmri_tle(atl)

    [confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();
    
    allsub_mean_eeg = zeros(5,sum(1:atlases(atl).regions-1));
    allsub_mean_fMRI = zeros(1,sum(1:atlases(atl).regions-1));
    
    allgroup_mean_EEG = zeros(length(groups), 5,sum(1:atlases(atl).regions-1));
    allgroup_mean_fMRI = zeros(length(groups),sum(1:atlases(atl).regions-1));
    
    tot_length = 69;
    
    allsub_fMRI = zeros(tot_length, sum(1:atlases(atl).regions-1));
    allsub_EEG = zeros(tot_length, 5, sum(1:atlases(atl).regions-1));
    group_idx = zeros(tot_length, 1);
    conf_idx = zeros(tot_length, 1);
    count_sub = 0;
    datamatrix = zeros(2, 3, 5, sum(1:atlases(atl).regions-1));

    for co = 1:length(confs)
        group_eegfmri = cell(3,1);
        
        for j = 1:length(groups)
            
            if(co==1)
                load([serialized_path 'eeg-fmri_' confs_label{co} '_' groups{j} '_connectomes_' atlases(atl).name '_scrubbed_truncTo5min'])
            else
                load([serialized_path 'eeg-fmri_' confs_label{co} '_' groups{j} '_connectomes_' atlases(atl).name '_scrubbed'])
            end

            conn_vec_length = length(subj(1).sess(1).fMRI);

            sub_mean_eeg = zeros(5,conn_vec_length);
            sub_mean_fMRI = zeros(1,conn_vec_length);
            
            subjects_eegfmri = zeros(length(subj),5);
            
            for i=1:length(subj)
                %increment total subject counter;
                count_sub = count_sub +1;
                group_idx(count_sub) = j;
                conf_idx(count_sub) = co;

                sess_mean_eeg = zeros(5,conn_vec_length);
                sess_mean_fMRI = zeros(1,conn_vec_length);

                %get session mean (average static of each session)
                for s = 1:length(subj(i).sess) 
                    sess_mean_fMRI = sess_mean_fMRI+subj(i).sess(s).fMRI'/length(subj(i).sess);
                    
                    for conn_type_it = 1:1 %cohi only length(subj(i).sess(s).EEG)
                        for b = 1:length(bands)
                            sess_mean_eeg(b, :) = sess_mean_eeg(b, :)+subj(i).sess(s).EEG(conn_type_it).bands(b).conn'/length(subj(i).sess);
                        end
                    end
                end

                %get subject mean per band
                for b = 1:length(bands)
                    sub_mean_eeg(b, :) = sub_mean_eeg(b, :) + sess_mean_eeg(b, :)/length(subj);
                    
                    allgroup_mean_EEG(j, b, :) = squeeze(allgroup_mean_EEG(j, b, :))' + sess_mean_eeg(b, :);
                    allsub_mean_eeg(b, :) = allsub_mean_eeg(b, :)+sess_mean_eeg(b, :);
                    allsub_EEG(count_sub, b, :) = sess_mean_eeg(b, :);
                    
                    subjects_eegfmri(i, b) = corr(squeeze(sess_mean_eeg(b, :))', sess_mean_fMRI');
                end
                
                allsub_fMRI(count_sub,:) = sess_mean_fMRI;
                
                allgroup_mean_fMRI(j, :) = allgroup_mean_fMRI(j, :) + sess_mean_fMRI;
                sub_mean_fMRI = sub_mean_fMRI + sess_mean_fMRI/length(subj);
                allsub_mean_fMRI = allsub_mean_fMRI+sess_mean_fMRI;
            end
            
            group_eegfmri{j} = subjects_eegfmri;

            for b = 1:length(bands)
                datamatrix(co, j, b) = corr(sub_mean_fMRI', sub_mean_eeg(b,:)');
            end
        end
    end
    
    %average group labels over all configurations
    for j = 1:length(groups)
        allgroup_mean_fMRI(j,:) = allgroup_mean_fMRI(j,:)/sum(group_idx==j);
        allgroup_mean_EEG(j,:,:) = allgroup_mean_EEG(j,:,:)/sum(group_idx==j);
    end
   
end