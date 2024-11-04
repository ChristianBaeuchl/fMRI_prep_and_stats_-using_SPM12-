%% Skull stripping
% performs skull stripping of anatomical images
% requires a bias-corrected anatomical image and 3 native-space tissue
% segmented images (GM, WM and CSF) of every subject
%
% by Christian Baeuchl - 18.02.2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear variables, close windows and clear the command line
clear , clc
tStart = tic;

%% Full path of experiment folder
path =  'path to your preprocessed data';
IDs  = cellstr(num2str(load('List.txt')));

%% skull stripping threshold
threshold = 0.5;

%% Initialise SPM defaults
spm('defaults','fmri');
spm_jobman('initcfg');

%-----------------------------------------------------------------------
% loop over IDs
for i = 1:numel(IDs)
    
    % load functional data for each session
    bias_corr_img =  spm_select('FPList', fullfile(path,[IDs{i} '_1'],'T1'), '^ms.*\.nii$');
    GM_seg_img    =  spm_select('FPList', fullfile(path,[IDs{i} '_1'],'T1'), '^c1s.*\.nii$');
    WM_seg_img    =  spm_select('FPList', fullfile(path,[IDs{i} '_1'],'T1'), '^c2s.*\.nii$');
    CSF_seg_img   =  spm_select('FPList', fullfile(path,[IDs{i} '_1'],'T1'), '^c3s.*\.nii$');
    
    input_imgs = {bias_corr_img;GM_seg_img;WM_seg_img;CSF_seg_img};
    
    matlabbatch{1}.spm.util.imcalc.input = input_imgs;
    matlabbatch{1}.spm.util.imcalc.output = [IDs{i} '_T1_skull_stripped.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr(fullfile(path,[IDs{i} '_1'],'T1'));
    matlabbatch{1}.spm.util.imcalc.expression = ['i1.*((i2+i3+i4)>' num2str(threshold) ')'];
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    %% SAVE AND RUN JOB
    %------------------
    save(fullfile(path,[IDs{i} '_1'],'T1',[IDs{i} '_skullstrip_batch.mat']),'matlabbatch');
    try
        spm_jobman('run', matlabbatch, '' );
        fprintf('Skull stripping performed on subjects %s anatomical T1 image\n',IDs{i})
    catch
        err_msg{i} = sprintf(['Could not run stats for subject/session: ' IDs{i}]);
    end
    
    % clean up variables
    clear bias_corr_img GM_seg_img WM_seg_img CSF_seg_img input_imgs
    
end
%%

tEnd = toc(tStart);
fprintf('\n ... done!\n');
fprintf('\nProcess finished in %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

