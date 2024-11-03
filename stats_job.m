%% stats_C1_job
% performs first level statistics
%
% Christian Baeuchl | last edited 14.05.2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function stats_job(datapath,statspath,IDs,SPECs)

% initialize matlabbatch cells
matlabbatch  = cell(1,4);
batch_matrix = cell(numel(IDs),4);

% loop over IDs
for i = 1:numel(IDs)

    % load functional data for each session
    funct_data =  spm_select('FPList',fullfile(datapath,num2str(IDs(i)),SPECs.functfold), SPECs.functpref);

    type = exist(fullfile(statspath,num2str(IDs(i)), SPECs.statsdir),'dir');
    if type ~= 7
        mkdir(fullfile(statspath,num2str(IDs(i))),SPECs.statsdir);
    end

    reg_file   = spm_select('FPList',fullfile(datapath,num2str(IDs(i)),SPECs.functfold), SPECs.mp);
    multi_cond = fullfile(SPECs.mcpath,[num2str(IDs(i)) SPECs.mc]);
    mc         = load(multi_cond);


    %% Display which subject and sequence is being processed
    fprintf('Starting 1st Level statistics for subject "%d" and the "%s" task of the "%s" session\n',...
        IDs(i),SPECs.task,SPECs.sess);

    %% CHANGE WORKING DIRECTORY (useful for .ps only)
    %-----------------------------------------------
    matlabbatch{1}.cfg_basicio.cfg_cd.dir = cellstr(fullfile(statspath,num2str(IDs(i))));

    %% FMRI model specification
    matlabbatch{2}.spm.stats.fmri_spec.dir            = cellstr(SPECs.statsdir);
    matlabbatch{2}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{2}.spm.stats.fmri_spec.timing.RT      = SPECs.rt;
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    matlabbatch{2}.spm.stats.fmri_spec.sess.scans     = cellstr(funct_data);
    matlabbatch{2}.spm.stats.fmri_spec.sess.cond      = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{2}.spm.stats.fmri_spec.sess.multi     = cellstr(multi_cond);
    matlabbatch{2}.spm.stats.fmri_spec.sess.regress   = struct('name', {}, 'val', {});
    matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = cellstr(reg_file);
    matlabbatch{2}.spm.stats.fmri_spec.sess.hpf       = 128;
    matlabbatch{2}.spm.stats.fmri_spec.fact           = struct('name', {}, 'levels', {});
    matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{2}.spm.stats.fmri_spec.volt           = 1;
    matlabbatch{2}.spm.stats.fmri_spec.global         = 'None';
    matlabbatch{2}.spm.stats.fmri_spec.mthresh        = SPECs.thresh;
    matlabbatch{2}.spm.stats.fmri_spec.mask           = {SPECs.brainmask};
    matlabbatch{2}.spm.stats.fmri_spec.cvi            = 'FAST';

    %% Model estimation
    matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File',...
        substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.fmri_est.write_residuals  = 0;
    matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;

    %% Contrast manager
    matlabbatch{4}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File',...
        substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

    %%
    if strcmp(SPECs.task,'GNG')
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.name    = 'Go (correct) > Nogo (correct)';
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0 0];
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{2}.tcon.name    = 'Nogo (correct) > Go (correct)';
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0 0];
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{3}.tcon.name    = 'Go (correct) | main effect';
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [1 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{4}.tcon.name    = 'Nogo (correct) | main effect';
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 0];
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

        if ~ismember(IDs(i),SPECs.excludeSubs)
            matlabbatch{4}.spm.stats.con.consess{5}.tcon.name    = 'Nogo (error) | main effect';
            matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights = [0 0 1 0];
            matlabbatch{4}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        end
        matlabbatch{4}.spm.stats.con.delete = 0;

        %%
    elseif strcmp(SPECs.task,'SST')
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.name    = 'Go (correct) > Stop (correct)';
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0 0];
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{2}.tcon.name    = 'Stop (correct) > Go (correct)';
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0 0];
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{3}.tcon.name    = 'Go (correct) | main effect';
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [1 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{4}.tcon.name    = 'Stop (correct) | main effect';
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 0];
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

        if ~ismember(IDs(i),SPECs.excludeSubs)
            matlabbatch{4}.spm.stats.con.consess{5}.tcon.name    = 'Stop (correct) > Stop (error)';
            matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights = [0 1 -1 0];
            matlabbatch{4}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

            matlabbatch{4}.spm.stats.con.consess{6}.tcon.name    = 'Stop (error) > Stop (correct)';
            matlabbatch{4}.spm.stats.con.consess{6}.tcon.weights = [0 -1 1 0];
            matlabbatch{4}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

            matlabbatch{4}.spm.stats.con.consess{7}.tcon.name    = 'Stop (error) | main effect';
            matlabbatch{4}.spm.stats.con.consess{7}.tcon.weights = [0 0 1 0];
            matlabbatch{4}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
        end
        matlabbatch{4}.spm.stats.con.delete = 0;

        %%
    elseif strcmp(SPECs.task,'STR')
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.name    = 'Congruent (correct) > Incongruent (correct)';
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0 0];
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{2}.tcon.name    = 'Incongruent (correct) > Congruent (correct)';
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0 0];
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{3}.tcon.name    = 'Congruent (correct) | main effect';
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [1 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

        matlabbatch{4}.spm.stats.con.consess{4}.tcon.name    = 'Incongruent (correct) | main effect';
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 0];
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        %
        n = 5;
        weights = zeros(1,4);
        indxCE = ismember(mc.names,'Congruent errors');
        indxIE = ismember(mc.names,'Incongruent errors');
        if any(indxCE) && any(indxIE)
            weights(find(indxCE)) = 1;
            weights(find(indxIE)) = -1;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'Congruent errors > Incongruent errors';
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            n = n+1;
            weights = -weights;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'Icongruent errors > Congruent errors';
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            n = n+1;
            weights = zeros(1,4);
        end
        if any(indxIE)
            weights(find(indxIE)) = 1;
            weights(2) = -1;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'Incongruent error > Icongruent correct';
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            n = n+1;
            weights = zeros(1,4);
            weights(find(indxIE)) = 1;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'Incongruent errors | main effect';
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            n = n+1;
            weights = zeros(1,4);
        end
        if any(indxCE)
            weights(find(indxCE)) = 1;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'Congruent errors | main effect';
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
            matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            n = n+1;
        end
        matlabbatch{4}.spm.stats.con.delete = 0;
        %
        if strcmp(SPECs.taskvar,'GLM02')
            weights = zeros(1,6);
            indxPIC = ismember(mc.names,'Post-incongruent correct trials');
            indxPIE = ismember(mc.names,'Post-incongruent errors');
            if any(indxPIE) && any(indxPIC)
                weights(find(indxPIE)) = 1;
                weights(find(indxPIC)) = -1;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'Post-icongruent (error) > Post-incongruent (correct)';
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n+1;
                weights = -weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'Post-incongruent (correct) > Post-icongruent (error)';
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            end
        end

        %%
    elseif strcmp(SPECs.task,'CDM')
        if strcmp(SPECs.taskvar,'GLM01') || strcmp(SPECs.taskvar,'GLM02')
            matlabbatch{4}.spm.stats.con.consess{1}.tcon.name    = 'Stimulus Decision';
            matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
            matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        end
        if strcmp(SPECs.taskvar,'GLM01')
            matlabbatch{4}.spm.stats.con.consess{2}.tcon.name    = 'Stimulus Decision | PM  Decision Value';
            matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
            matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

            matlabbatch{4}.spm.stats.con.consess{3}.tcon.name    = 'Stimulus Decision | PM  Decision Value | negative contrast';
            matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [0 -1 0];
            matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            matlabbatch{4}.spm.stats.con.delete = 0;
        elseif strcmp(SPECs.taskvar,'GLM02')
            matlabbatch{4}.spm.stats.con.consess{2}.tcon.name    = 'Stimulus Decision | PM  Shortterm Consequences';
            matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
            matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

            matlabbatch{4}.spm.stats.con.consess{3}.tcon.name    = 'Stimulus Decision | PM  Longterm Consequences';
            matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
            matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            matlabbatch{4}.spm.stats.con.delete = 0;
        end
        if strcmp(SPECs.taskvar,'GLM03')
            n = 1;
            weights = zeros(1,8);
            indxSCSL = ismember(mc.names,'Self Control Success | Longterm');
            indxSCSS = ismember(mc.names,'Self Control Success | Shortterm');
            if any(indxSCSL) && any(indxSCSS)
                weights(find(indxSCSL)) = 1;
                weights(find(indxSCSS)) = -1;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'SC Success LongTerm > SC Success ShortTerm ';
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n+1;
                weights = -weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'SC Success ShortTerm > SC Success LongTerm ';
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n+1;
                weights = zeros(1,8);
            end
            indxSCFL = ismember(mc.names,'Self Control Failure | Longterm');
            indxSCFS = ismember(mc.names,'Self Control Failure | Shortterm');
            if any(indxSCFL) && any(indxSCFS)
                weights(find(indxSCFL)) = 1;
                weights(find(indxSCFS)) = -1;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'SC Failure LongTerm > SC Failure ShortTerm ';
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n+1;
                weights = -weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.name    = 'SC Failure ShortTerm > SC Failure LongTerm ';
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.weights = weights;
                matlabbatch{4}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            end
            matlabbatch{4}.spm.stats.con.delete = 0;
        end
    end

    %% CREATE MATLABBATCH MATRIX AND SAVE JOB
    batch_matrix(i,:) = matlabbatch;
    save(fullfile(statspath,num2str(IDs(i)),[num2str(IDs(i)) '_' SPECs.statsdir '.mat']),'matlabbatch');
    matlabbatch       = cell(1,4);

end

%% RUN JOB [using parallelisation]
%------------------
parfor i = 1:numel(IDs)
    try
        spm_jobman('serial', batch_matrix(i,:), '' );
    catch
        err_msg{i} = sprintf(['Could not run stats for subject/session: ' num2str(IDs(i))]);
    end
end

%% DEBLANK ERROR MESSAGE VARIABLE
if ~exist('err_msg','var')
    fprintf('All datasets were analyzed\n')
else
    err_msg = err_msg(~cellfun('isempty',err_msg));
    fprintf('\n------------------------------------------------------------------------------\n');
    fprintf(2,'!!!At least one data set was not preprocessed: check out variable "err_msg"!!!\n');
    fprintf('------------------------------------------------------------------------------\n');
    assignin('base','err_msg',err_msg);
end
