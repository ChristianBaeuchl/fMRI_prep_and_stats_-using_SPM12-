%% DARTEL normalization to MNI space
%
% Normalizes (writes) the functional EPIs to MNI space
% Requires that a group template (or templates) as well as subject-specific
% flow fields already exist. To run sucessfully the path structure that the
% code requires must be adhered to (or these parts need to be changed
% within the code)
%
% Calls upon the normalize2mni_job_200113 file to run the procedure and
% uses MATLABs parallelization toolbox to speed-up the job.
%
% by Christian Baeuchl | last edited on 11.12.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear, clc    % clear workspace and command window 
tStart = tic; % start counting the code running time

%%
% voxel size of EPI images
VARS.vox_size = [3 3 3];
% smoothing vector
VARS.smooth_vec = [8 8 8];
% bounding box
VARS.b_box = [-78 -112 -50; 78 76 85];
% functional data folder
VARS.exp           = 'task name';
VARS.se1           = 'session 1 name';
VARS.se2           = 'session 2 name';
VARS.dat           = 'date of template creation';
VARS.epi_folder    = 'task folder name';
VARS.struct_folder = 'structural folder name';

%% Full path of experiment folder
datapath   = ['path to your preprocessed data' VARS.se1 '- uses session information'];
scriptpath = 'path to your DARTEL_normalize2mni_job code file';

%%  Subject folders
IDs = cellstr(num2str(importdata(['./Participant_lists/List_' VARS.exp '_' VARS.se1  '.txt'])));

%% Initialise SPM defaults
spm('Defaults','fMRI');
spm_jobman('initcfg');

%% Full path of YA template
template_name   = ['Template_' VARS.exp '_' VARS.se2 '_' VARS.dat '_6.nii'];
template_folder = IDs{1};
template_path   = [datapath filesep template_folder filesep VARS.struct_folder filesep template_name];
flow_name       = ['^u_.*\Template_' VARS.exp '_' VARS.se2 '_' VARS.dat '.nii$'];

normalize2mni_job_C1_231211(datapath,scriptpath,IDs,template_path,flow_name,VARS); % now with parallelization

%%
tEnd = toc(tStart);
fprintf('\n ... done!\n');
fprintf('\nProcess finished in %d hours, %d minutes, %d seconds and %d milliseconds\n',...
    floor(tEnd/(60*60)),floor(rem((tEnd/60),60)),floor(rem(tEnd,60)),round((tEnd-floor(tEnd))*1000));

