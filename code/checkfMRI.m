% Get permuted split-half of fMRI connectomes
% SI Table 3 manuscript
%
% 2022-08-22 Jonathan Wirsich
[confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

perm_iters = 5000;
perm_iters_conf = zeros(length(confs), 1);

for atl=1:1
    
    conn_group = cell(length(confs),1);
    conn_group_splithalf = cell(length(confs),2,perm_iters);
    
    for j = 1:length(groups)
        disp(groups{j})
        for co = 1:length(confs)    
            %load struct subj
            if(co==1) %Geneva dataset load cut to 5 minutes
                load([serialized_path 'eeg-fmri_' confs{co} '_'  groups{j} '_connectomes_' atlases(atl).name '_scrubbed_truncTo5min'])
            else
                load([serialized_path 'eeg-fmri' '_' confs{co} '_' groups{j} '_connectomes_' atlases(atl).name '_scrubbed'])
            end

            conn_vec_length = length(subj(1).sess(1).fMRI);
            sub_mean = zeros(1,conn_vec_length);

            %for n<16 fully sample all combinations with nchoosek(n, round(n/2));
            if length(subj) < 16
                nchoosek_subj = nchoosek(1:length(subj), round(length(subj)/2));
                %get only half of the combination as rest of the combination
                %will be asigned to second group (limitation: unbalanced groups
                %in case of unpair group size)
                perm_iters = floor(length(nchoosek_subj)/2);
            end
            %save determined iterations for display
            perm_iters_conf(co) = perm_iters;

            for iter = 1:perm_iters
                sub_mean_split = zeros(2, conn_vec_length);
                %for n<16 fully sample all combinations with nchoosek(n, round(n/2));
                if length(subj) < 16
                    %get permutation and rest of the set
                    perm_subidx = [nchoosek_subj(iter, :) setdiff(1:length(subj), nchoosek_subj(iter, :))];                
                else
                    perm_subidx = randperm(length(subj));
                end

                %count loop iteration to know when to split half
                count = 0;
                for i = perm_subidx
                    sess_mean = zeros(1,conn_vec_length);
                    for s = 1:length(subj(i).sess)
                        sess_mean = sess_mean+subj(i).sess(s).fMRI'/length(subj(i).sess);
                    end

                    %needs to be calculated only one time
                    if iter == 1
                        sub_mean = sub_mean + sess_mean/length(subj);
                    end

                    count = count +1;
                    %split at count-1 (group 1 has more subjects in cae of
                    %unpair splits / rounding correct in case of pair splits)
                    split_idx = round((count-1)/length(subj));
                    if split_idx == 0
                        sub_mean_split(split_idx+1,:) = sub_mean_split(split_idx+1,:) ...
                            + sess_mean/(round(length(subj)/2));
                    else
                        sub_mean_split(split_idx+1,:) = sub_mean_split(split_idx+1,:) ...
                            + sess_mean/(floor(length(subj)/2));
                    end
                end

                conn_group_splithalf{co, 1, iter} = model.Connectome(sub_mean_split(1,:),  atlases(atl).regions);
                conn_group_splithalf{co, 2, iter} = model.Connectome(sub_mean_split(2,:),  atlases(atl).regions);
            end
            conn_group{co} = model.Connectome(sub_mean,  atlases(atl).regions);
        end

        disp('Inter split')
        for c1 = 1:length(confs)
            perm_iters = perm_iters_conf(c1);
            corr_perm = zeros(1, perm_iters);
            for iter = 1:perm_iters
                con1 = conn_group_splithalf{c1, 1, iter}.vec';
                con2 = conn_group_splithalf{c1, 2, iter}.vec';
                corr_perm(iter) = corr(con1, con2);
            end
            disp([confs_label{c1} ' ' num2str(mean(corr_perm))]);
        end

        disp('Inter dataset')
        for c1 = 1:length(confs)-1
            for c2 = c1+1:length(confs)
                con1 = conn_group{c1}.vec';
                con2 = conn_group{c2}.vec';
                disp([confs_label{c1} '-' confs_label{c2} ' ' num2str(corr(con1, con2))]);
            end
        end
    end
    
end
