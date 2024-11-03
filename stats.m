%% 1st level fMRI stats (SPM12)
% performs first level statistics of the 4 fMRI tasks:
%
%   Complex Decision Making (CDM)
%   Go/Nogo task (GNG)
%   Stop-Signal task (SST)
%   Counting Stroop task (STR)
%
% last edited by Christian Baeuchl | 22.05.2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear , clc, close
tStart = tic;

%% define specifications for statistical analysis
SPECs.task      = 'GNG';
SPECs.taskvar   = 'GLM01';
SPECs.glmsuffix = '_mk2021';
SPECs.idsuffix  = '';
SPECs.sess      = '01_BL';
SPECs.functfold = 'GoNogo_EPI';
SPECs.functpref = '^swa.*\.nii$';
SPECs.statsdir  = ['stats_' SPECs.task '_' SPECs.taskvar];
SPECs.stats     = 'SPM.mat';
SPECs.mp        = '^rp.*\.txt$';
SPECs.rt        = 2;
SPECs.thresh    = 0.8;
SPECs.mc        = ['_' SPECs.task '_' SPECs.sess(4:end) '.mat'];
SPECs.mcpath    = 'path to your "multiple conditions" a.k.a. "stimulus onset" file';

%% On which machine does the analysis run on?
machine = 'office_pc'; %'hexacore'; %'office_pc'; %
%%
if strcmp(machine,'office_pc')
    SPECs.brainmask = 'P:\211\projektmitarbeiter_intern\Daten\06_MRI_analyses\06_Analyses_material\Brainmasks\binary_MNI152_T1_3mm_brain.nii';
elseif strcmp(machine,'hexacore')
    SPECs.brainmask = 'D:\Programs\MATLAB\Toolboxes\spm12\spm12\toolbox\FieldMap\brainmask.nii';
end

%% Data and stats path
datapath  = ['path to your preprocessed data' SPECs.sess '- includes session information'];
statspath = ['path to where your stats data should be stored' SPECs.sess '- includes session information'];

IDs = sort(importdata(['.\Participant_lists\List_' SPECs.task '_' SPECs.sess SPECs.idsuffix '.txt']));

%% Initialise SPM defaults
spm('defaults','fmri');
spm_jobman('initcfg');

%% Call the job function
stats_job(datapath,statspath,IDs,SPECs);

%%
tEnd = toc(tStart);
fprintf('\n ... done!\n');
fprintf('\nProcess finished in %d hours, %d minutes, %d seconds and %d milliseconds\n',...
    floor(tEnd/(60*60)),floor(rem((tEnd/60),60)),floor(rem(tEnd,60)),round((tEnd-floor(tEnd))*1000));

