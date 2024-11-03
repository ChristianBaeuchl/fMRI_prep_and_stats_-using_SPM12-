%% DARTEL - create group template
%
% Christian Baeuchl | last edited 07.12.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc    % clear workspace and command window 
tStart = tic; % start counting the code running time

%% Full path of experiment folder
session    = {'list session names here'};
datapath   = ['path to your preprocessed data' session{2} ' includes session information'];
scriptpath = 'path to your DARTEL_create_template_job script';
IDs.structfolder  = 'name of folder with your structdural image';
IDs.templatename  = ['Template_' session{3} '_' session{1}];
IDs.skullstripped = 0; % if '1' % uses the skull-stripped structural image if '1', '0': otherwise

%%  IDs
IDs.list = num2str(importdata(['./Participant_lists/List_' session{3} '_' session{2} '.txt']));

%% Initialise SPM defaults
spm('Defaults','fMRI');
spm_jobman('initcfg');

DARTEL_create_template_job(session,datapath,scriptpath,IDs);

%%
tEnd = toc(tStart);
fprintf('\n ... done!\n');
fprintf('\nProcess finished in %d hour(s), %d minute(s), %d second(s) and %d milliseconds\n',...
    floor(tEnd/(60*60)),floor(rem((tEnd/60),60)),floor(rem(tEnd,60)),round((tEnd-floor(tEnd))*1000));
