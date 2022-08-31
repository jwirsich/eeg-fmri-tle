% init global variables for EEG-fMRI connectomes
% set serialized_path to cloned repo from https://github.com/jwirsich/eeg-fmri-tle/data
%
% 2022-08-22 Jonathan Wirsich
function [confs, confs_label, groups, eeg_conn_types, bands, atlases, serialized_path] = init_serialized_conf_tle()

    co = 0;
    co = co +1;
    confs{co} = '256Ch3T';
    confs_label{co} = '256Ch3T';
    co = co +1;
    confs{co} = '64Ch3T';
    confs_label{co} = '64Ch3T';

    groups{1} = 'controls';
    groups{2} = 'ltle';
    groups{3} = 'rtle';
    
    %should probably be read out from matfile
    eeg_conn_types = cell(3, 1);
    eeg_conn_types{1} = 'cohi';

    bands{1} = 'delta';
    bands{2} = 'theta';
    bands{3} = 'alpha';
    bands{4} = 'beta';
    bands{5} = 'gamma';
    
    %consider making an object
    atlases = struct('name', [], 'regions', []);
    atlases(1).name = 'desikan';
    atlases(1).regions = 68;

    reflect_folder = fileparts(mfilename('fullpath'));
    
    serialized_path= [reflect_folder filesep '..' filesep '..' filesep 'data' filesep];

end