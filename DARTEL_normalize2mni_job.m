%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

function normalize2mni_job_C1_231211(datapath,scriptpath,IDs,template_path,flow_name,VARS)


%%
matlabbatch{1}.spm.tools.dartel.mni_norm.template = cellstr(template_path);

for i = 1:numel(IDs)
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(i).flowfield = cellstr(spm_select('FPList',[datapath filesep IDs{i} filesep VARS.struct_folder],flow_name));
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(i).images = cellstr(spm_select('FPList',fullfile(datapath,IDs{i},VARS.epi_folder),'^au.*\.nii$'));
    if isempty(matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(i).images)
        matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj(i).images = cellstr(spm_select('FPList',fullfile(datapath,IDs{i},VARS.epi_folder),'^ar.*\.nii$'));
    end
end

matlabbatch{1}.spm.tools.dartel.mni_norm.vox      = VARS.vox_size;
matlabbatch{1}.spm.tools.dartel.mni_norm.bb       = VARS.b_box;
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm     = VARS.smooth_vec;

%
save(fullfile(scriptpath,['Normalize2MNI_C1_' VARS.se2 '_' VARS.exp '_' datestr(now,'yymmdd') '.mat']),'matlabbatch');


%% split-up matlabbatch for parallelization

% determine how many worker are on current machine
myCluster = parcluster('local');
myCluster.NumWorkers;
fprintf('The number of workers for the normalization jobs is: %s\n', num2str(myCluster.NumWorkers));

if size(matlabbatch,2) == 1
    
    numsubs_01 = size(matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj,2);
    binwidth_01 = numsubs_01/myCluster.NumWorkers;
    
    if binwidth_01 == round(binwidth_01)
        temp = matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj;
        batchtemp = matlabbatch(1,1);
        batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        n = 1;
        for i = 1:myCluster.NumWorkers
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:binwidth_01*i);
            batchmatrix(i,:) = batchtemp;
            n = (binwidth_01*i) + 1;
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        end
        
    elseif binwidth_01 ~= round(binwidth_01)
        binwidth_01 = floor(binwidth_01);
        subsdiff_01 = numsubs_01 - (binwidth_01*myCluster.NumWorkers);
        temp = matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj;
        batchtemp = matlabbatch(1,1);
        batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        n = 1;
        for i = 1:myCluster.NumWorkers
            if i == myCluster.NumWorkers
                batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:((binwidth_01*i) + subsdiff_01));
                batchmatrix(i,:) = batchtemp;
                break
            end
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:binwidth_01*i);
            batchmatrix(i,:) = batchtemp;
            n = (binwidth_01*i) + 1;
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        end
    end
    
elseif size(matlabbatch,2) == 2
    
    numsubs_01 = size(matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj,2);
    binwidth_01 = numsubs_01/myCluster.NumWorkers;
    numsubs_02 = size(matlabbatch{1,2}.spm.tools.dartel.mni_norm.data.subj,2);
    binwidth_02 = numsubs_02/myCluster.NumWorkers;
    
    if binwidth_01 == round(binwidth_01)
        temp = matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj;
        batchtemp = matlabbatch(1,1);
        batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        n = 1;
        for i = 1:myCluster.NumWorkers
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:binwidth_01*i);
            batchmatrix_01(i,:) = batchtemp;
            n = (binwidth_01*i) + 1;
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        end
        
    elseif binwidth_01 ~= round(binwidth_01)
        binwidth_01 = floor(binwidth_01);
        subsdiff_01 = numsubs_01 - (binwidth_01*myCluster.NumWorkers);
        temp = matlabbatch{1,1}.spm.tools.dartel.mni_norm.data.subj;
        batchtemp = matlabbatch(1,1);
        batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        n = 1;
        for i = 1:myCluster.NumWorkers
            if i == myCluster.NumWorkers
                batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:((binwidth_01*i) + subsdiff_01));
                batchmatrix_01(i,:) = batchtemp;
                break
            end
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:binwidth_01*i);
            batchmatrix_01(i,:) = batchtemp;
            n = (binwidth_01*i) + 1;
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        end
    end
    
    if binwidth_02 == round(binwidth_02)
        temp = matlabbatch{1,2}.spm.tools.dartel.mni_norm.data.subj;
        batchtemp = matlabbatch(1,2);
        batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        n = 1;
        for i = 1:myCluster.NumWorkers
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:binwidth_02*i);
            batchmatrix_02(i,:) = batchtemp;
            n = (binwidth_02*i) + 1;
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        end
        
    elseif binwidth_02 ~= round(binwidth_02)
        binwidth_02 = floor(binwidth_02);
        subsdiff_02 = numsubs_02 - (binwidth_02*myCluster.NumWorkers);
        temp = matlabbatch{1,2}.spm.tools.dartel.mni_norm.data.subj;
        batchtemp = matlabbatch(1,2);
        batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        n = 1;
        for i = 1:myCluster.NumWorkers
            if i == myCluster.NumWorkers
                batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:((binwidth_02*i) + subsdiff_02));
                batchmatrix_02(i,:) = batchtemp;
                break
            end
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data.subj = temp(n:binwidth_02*i);
            batchmatrix_02(i,:) = batchtemp;
            n = (binwidth_02*i) + 1;
            batchtemp{1,1}.spm.tools.dartel.mni_norm.data = [];
        end
    end
    
    batchmatrix = [batchmatrix_01;batchmatrix_02];
end


%% RUN JOB
%------------------
parfor i = 1:size(batchmatrix,1)
    try
        spm_jobman('serial', batchmatrix(i,:), '' );
    catch
        err_msg{i} = sprintf(['Could not run normalization for cluster #: ' num2str(i)]);
    end
end

%% DEBLANK ERROR MESSAGE VARIABLE
if ~exist('err_msg','var')
    fprintf('All datasets were analyzed\n')
else
    err_msg = err_msg(~cellfun('isempty',err_msg));
    fprintf('\n--------------------------------------------------------------------------------------\n');
    fprintf(2,'!!!At least one cluster data set was not normalized: check out variable "err_msg"!!!\n');
    fprintf('----------------------------------------------------------------------------------------\n');
    assignin('base','err_msg',err_msg);
end

end
