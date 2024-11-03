%% preprocessing_job script
% The function performs spatial preprocessing for each subject and assigns
% appropriate values and files. The spm batch jobs are created and run *PARALLELLY*.
%
% Parallel processes run in the PARFOR LOOP
%
% Christian Baeuchl | latest update 30.10.2023

function preprocessing_job(datapath,IDs,VDM,RUN,STC,IMG,TMP)

% get index for matlabbatch
switch STC.opt
    case 0
        if strcmp(IMG.Opt,'CoregSeg')
            num = [1 2 3 NaN 4 5];
        elseif strcmp(IMG.Opt,'NoCoregSeg')
            num = [1 2 3 4];
        elseif strcmp(IMG.Opt,'OnlyCoregSeg')
            num = [1 NaN NaN NaN 2 3];
        end
    case 1
        if strcmp(IMG.Opt,'CoregSeg')
            num = [1 2 4 3 5 6];
        elseif strcmp(IMG.Opt,'NoCoregSeg')
            num = [1 2 3 4 NaN NaN];
        end
    case 2
        if strcmp(IMG.Opt,'CoregSeg')
            num = [1 2 3 4 5 6];
        elseif strcmp(IMG.Opt,'NoCoregSeg')
            num = [1 2 3 4 NaN NaN];
        end
end

% initialize matlabbatch cells
matlabbatch = cell(1,nnz(~isnan(num)));
batch_matrix = cell(numel(IDs),nnz(~isnan(num)));

% loop over IDs
for i = 1:numel(IDs)

    % load data
    funct_data = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.EPI),['^' num2str(IDs(i)) '.*\_4D.nii$']));
    structural = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.T1), '^s.*\.nii$'));
    if isempty(structural{1})
        structural = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.T1),'^MF\d.*.nii'));
    end

    if strcmp(IMG.Opt,'OnlyCoregSeg')
        meanimg = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.EPI),'^mean.*\_4D.nii$'));
    end

    % load fieldmap images
    phase     = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.Ph),'^s\d.*.nii'));
    if isempty(phase{1})
        phase = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.Ph),'^MF\d.*.nii'));
    end
    magnitude     = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.Mag),'^s\d.*01.nii'));
    if isempty(magnitude{1})
        magnitude = cellstr(spm_select('FPList', fullfile(datapath,num2str(IDs(i)),IMG.Mag),'^MF\d.*01.nii'));
    end

    %% display which subject and sequence is being processed
    fprintf('Preprocessing subject  - "%d" - %s | %s\n',IDs(i),IMG.NAM,IMG.Se1);

    %% CHANGE WORKING DIRECTORY (useful for .ps only)
    %-----------------------------------------------
    matlabbatch{num(1)}.cfg_basicio.cfg_cd.dir = cellstr(fullfile(datapath,num2str(IDs(i))));

    %% CALCULATE VDM
    if ~isnan(num(2))
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.phase      = phase;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.magnitude  = magnitude;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.et        = [VDM.short_te VDM.long_te];
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.maskbrain = 1;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir   = VDM.blip;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.tert      = VDM.total_rot;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.epifm     = 0;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.ajm       = 0;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.method = 'Mark3D';
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.fwhm   = 10;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.pad    = 0;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.ws     = 1;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.template = TMP.fieldmap;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.fwhm    = 5;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.nerode  = 2;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.ndilate = 4;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.thresh  = 0.5;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.reg     = 0.02;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).session.epi = funct_data;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).matchvdm    = 1;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).sessname    = 'session';
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).writeunwarped = 0;
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).anat        = '';
        matlabbatch{num(2)}.spm.tools.fieldmap.calculatevdm.subj(1).matchanat   = 0;
    end


    %% REALIGN & UNWARP
    if ~isnan(num(3))
        if STC.opt == 2 || STC.opt == 0
            matlabbatch{num(3)}.spm.spatial.realignunwarp.data(1).scans     = funct_data;
            matlabbatch{num(3)}.spm.spatial.realignunwarp.data(1).pmscan(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 1, Session 1)',...
                substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{1}));
        elseif STC.opt == 1
            matlabbatch{num(3)}.spm.spatial.realignunwarp.data(1).scans(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)',...
                substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
            matlabbatch{num(3)}.spm.spatial.realignunwarp.data(1).pmscan(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 1, Session 1)',...
                substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{1}));
        end
        matlabbatch{num(3)}.spm.spatial.realignunwarp.eoptions.quality = 1;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.eoptions.sep     = 2;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.eoptions.fwhm    = 5;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.eoptions.rtm     = 0;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.eoptions.einterp = 7;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.eoptions.ewrap   = RUN.ped;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.eoptions.weight  = '';
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.jm  = 0;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uwroptions.uwwhich  = [2 1];
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uwroptions.rinterp  = 7;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uwroptions.wrap     = RUN.ped;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uwroptions.mask     = 1;
        matlabbatch{num(3)}.spm.spatial.realignunwarp.uwroptions.prefix   = 'u';
    end

    %% Slice Time Correction
    if ~isnan(num(4))
        if STC.opt ~= 0
            if STC.opt == 1
                matlabbatch{num(4)}.spm.temporal.st.scans{1} = funct_data;
            elseif STC.opt == 2
                matlabbatch{num(4)}.spm.temporal.st.scans{1}(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)',...
                    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
            end
            matlabbatch{num(4)}.spm.temporal.st.nslices  = STC.nslices;
            matlabbatch{num(4)}.spm.temporal.st.tr       = STC.tr;
            matlabbatch{num(4)}.spm.temporal.st.ta       = STC.ta;
            matlabbatch{num(4)}.spm.temporal.st.so       = STC.so;
            matlabbatch{num(4)}.spm.temporal.st.refslice = STC.refslice;
            matlabbatch{num(4)}.spm.temporal.st.prefix   = 'a';
        end
    end

    %% COREGISTRATION
    if strcmp(IMG.Opt,'OnlyCoregSeg')
        matlabbatch{num(5)}.spm.spatial.coreg.estimate.ref(1) = meanimg;
    else
        matlabbatch{num(5)}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image',...
            substruct('.','val', '{}',{num(3)}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
    end
    matlabbatch{num(5)}.spm.spatial.coreg.estimate.source = structural;
    matlabbatch{num(5)}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{num(5)}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{num(5)}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{num(5)}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{num(5)}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    %% SEGMENTATION
    matlabbatch{num(6)}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images',...
        substruct('.','val', '{}',{num(5)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{num(6)}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{num(6)}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{num(6)}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(1).tpm = TMP.tpm1;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(1).native = [1 1];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(2).tpm = TMP.tpm2;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(2).native = [1 1];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(3).tpm = TMP.tpm3;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(3).native = [1 1];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(4).tpm = TMP.tpm4;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(5).tpm = TMP.tpm5;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(6).tpm = TMP.tpm6;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{num(6)}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{num(6)}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{num(6)}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{num(6)}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{num(6)}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{num(6)}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{num(6)}.spm.spatial.preproc.warp.write = [0 1];


    %% CREATE MATLABBATCH MATRIX AND SAVE JOB
    batch_matrix(i,:) = matlabbatch;
    save(fullfile(datapath,num2str(IDs(i)),strcat(num2str(IDs(i)),'_preprocessing_',IMG.NAM,'_',datestr(now,'yymmdd'),'.mat')),'matlabbatch');
    matlabbatch = cell(1,nnz(~isnan(num)));

end

%% RUN JOB
%------------------
parfor i = 1:numel(IDs)
    try
        spm_jobman('serial', batch_matrix(i,:), '' );
    catch
        err_msg{i} = sprintf(['Could not run preprocessing for subject: ' num2str(IDs(i))]);
    end
end

%% DEBLANK ERROR MESSAGE VARIABLE
if ~exist('err_msg','var')
    fprintf('All datasets were preprocessed\n')
else
    err_msg = err_msg(~cellfun('isempty',err_msg));
    fprintf('\n------------------------------------------------------------------------------\n');
    fprintf(2,'!!!At least one data set was not preprocessed: check out variable "err_msg"!!!\n');
    fprintf('------------------------------------------------------------------------------\n');
    assignin('base','err_msg',err_msg);
end

