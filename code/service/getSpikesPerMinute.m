% load marked spikes from data
%
% conf: name of configuration ('eeg-fmri-geneva' or 'eeg-fmri-marseille')
% group: name of group ('lTLE' or 'rTLE')
% isVerbose: boolean print Spikes per minute (true/false)
%
% 2022-08-22 Jonathan Wirsich
function spkPerMin = getSpikesPerMinute(conf, group, isVerbose)

    [confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle();

    if strcmp(conf, '256Ch3T')
        c = 1;
    elseif strcmp(conf, '64Ch3T')
        c = 2;
    else
        error(['Configuration ' conf  ' not defined'])
    end
    
    if(isVerbose)
        disp('Verbose Spikes')
        disp(['conf: ' confs_label{c}])
    end

    if c == 1
        tmp = load([serialized_path 'eeg-fmri_' confs{c} '_tle_spikes.mat']);
        tr = 1.990;
        no_vol = 600;
    else
        tmp = load([serialized_path 'eeg-fmri_' confs{c} '_tle_spikes.mat']);
        tr = 3.6;
        no_vol = 345;
    end
    subj = tmp.subj;

    count_subj = 0;
    for it_subj = 1:length(subj)
        
       if strcmpi(subj(it_subj).group, group)
           count_subj = count_subj + 1;
           spikecount(count_subj) = length(subj(it_subj).spike_event);
           spikecountPerMin(count_subj) = length(subj(it_subj).spike_event)/(tr*no_vol/60);
           spikecount_5min(count_subj) = sum(subj(it_subj).spike_event<=150);
           spikecount_5minPerMin(count_subj) = sum(subj(it_subj).spike_event<=150)/(tr*150/60);

           if(isVerbose)
               if c == 1
                  disp([subj(count_subj).name  ': ' num2str(spikecountPerMin(count_subj))]) 
               else
                  disp([subj(count_subj).name  ': ' num2str(spikecount_5minPerMin(count_subj))])
               end
           end
       end
    end
    
    if c == 1
       spkPerMin = spikecount_5minPerMin';
    else
       spkPerMin = spikecountPerMin';
    end

end
