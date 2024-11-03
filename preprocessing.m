%% preprocessing script
%  This script defines paths, folders and parameters as input for the
%  preprocessing_job function (which is called upon at the end of this
%  script)
%
% Christian Baeuchl | latest update 30.10.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear,clc % clear workspace and command window

%% parameters for calculating VDM (voxel displacement map)
VDM.short_te  = 5.32;   % TE 1;
VDM.long_te   = 7.78;   % TE 2;
VDM.total_rot = 37.12;  % total readout time = # of echoes * echo spacing | trt = EPI factor * echo spacing
VDM.blip      = -1;     % Blip direction - phase encoding direction A >> P, blip = -1; P >> A, blip = 1;

%% parameters/options for slice time correction (STC)
STC.opt      = 2;  %  0 = no STC; 1 = STC before realignment; 2 = STC after realignment
STC.nslices  = 34; % number of slices 
STC.tr       = 2;  % TR = time of repetition
STC.ta       = STC.tr-(STC.tr/STC.nslices); % TA - time of acquisition
STC.so       = STC.nslices:-1:1; % slice acquision order (descending in this case)
STC.refslice = 17; % which slice should be the reference?

%% parameters/options for realign & unwarp
RUN.ped = [0 1 0]; % indicates phase encoding direction. [1 0 0] = X (L >> R), [0 1 0] = Y (A >> P); [0 0 1] = X & Y

%% Subjects and image folders (not full path)

%%  Image Folders
IMG.T1  = 'T1';                     % MP-RAGE
IMG.Mag = 'Fieldmap\Magnitude';     % Magnitude
IMG.Ph  = 'Fieldmap\Phase';         % Phase
IMG.EPI = 'Stroop_EPI';             % EPI Folder
IMG.NAM = 'STR';                    % Task abbreviation 
IMG.Se1 = 'Follow-up 3';            % Session information
IMG.Se2 = '02_FU3';                 % Session information
IMG.Opt = 'CoregSeg';               % 'OnlyCoregSeg'; 'NoCoregSeg';

%% Full path of experiment folder
datapath    = ['path to your functional data' IMG.Se2 '- contains session information'];
workstation = 'hexaCore'; % 'officePC'; % in case you run this code on different PCs with different paths to certain files (see down below)

%%  Subject IDs

%IDs = load(['./Participant_lists/List_' IMG.NAM '_' IMG.Se2 '.txt']);
IDs = [subject_01
       subject_02
       subject_03,...
       ];

%%  Template Paths
if strcmp(workstation,'officePC')
    TMP.fieldmap = {'path to your fieldmap template'};
    TMP.tpm1 =     {'path to your tissue probability map (TPM) #1: ...\TPM.nii,1'};
    TMP.tpm2 =     {'path to your tissue probability map (TPM) #2: ...\TPM.nii,2'};
    TMP.tpm3 =     {'path to your tissue probability map (TPM) #3: ...\TPM.nii,3'};
    TMP.tpm4 =     {'path to your tissue probability map (TPM) #4: ...\TPM.nii,4'};
    TMP.tpm5 =     {'path to your tissue probability map (TPM) #5: ...\TPM.nii,5'};
    TMP.tpm6 =     {'path to your tissue probability map (TPM) #6: ...\TPM.nii,6'};
elseif strcmp(workstation,'hexaCore')
    TMP.fieldmap = {'path to your fieldmap template'};
    TMP.tpm1 =     {'path to your tissue probability map (TPM) #1: ...\TPM.nii,1'};
    TMP.tpm2 =     {'path to your tissue probability map (TPM) #2: ...\TPM.nii,2'};
    TMP.tpm3 =     {'path to your tissue probability map (TPM) #3: ...\TPM.nii,3'};
    TMP.tpm4 =     {'path to your tissue probability map (TPM) #4: ...\TPM.nii,4'};
    TMP.tpm5 =     {'path to your tissue probability map (TPM) #5: ...\TPM.nii,5'};
    TMP.tpm6 =     {'path to your tissue probability map (TPM) #6: ...\TPM.nii,6'};
end
   
%% Initialise SPM defaults
spm('Defaults','fMRI');
spm_jobman('initcfg');

preprocessing_job(datapath,IDs,VDM,RUN,STC,IMG,TMP);
