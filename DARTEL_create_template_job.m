%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

function DARTEL_create_template_job(session,datapath,scriptpath,IDs)

%-----------------------------------------------------------------------
% Job saved on 06-Feb-2017 14:02:51 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.cfg_cd.dir  = cellstr(datapath);

%% Create Template Young Adults
for i = 1:size(IDs.list,1)
    GM_img      = cellstr(spm_select('FPList',[datapath filesep IDs.list(i,:) filesep IDs.structfolder], '^rc1.*\.nii$'));
    GM_img_name = cellstr(spm_select('List',  [datapath filesep IDs.list(i,:) filesep IDs.structfolder], '^rc1.*\.nii$'));
    if size(GM_img,1) > 1 && IDs.skullstripped
        GM_img = GM_img(endsWith(GM_img,'skull_stripped.nii'));
    elseif size(GM_img,1) > 1 && ~IDs.skullstripped
        GM_img = GM_img(startsWith(GM_img_name,'rc1s'));
    end
    matlabbatch{2}.spm.tools.dartel.warp.images{1}(i,1) = GM_img;
    
    WM_img      = cellstr(spm_select('FPList',[datapath filesep IDs.list(i,:) filesep IDs.structfolder], '^rc2.*\.nii$'));
    WM_img_name = cellstr(spm_select('List',  [datapath filesep IDs.list(i,:) filesep IDs.structfolder], '^rc2.*\.nii$'));
    if size(WM_img,1) > 1 && IDs.skullstripped
        WM_img = WM_img(endsWith(WM_img,'skull_stripped.nii'));
    elseif size(WM_img,1) > 1 && ~IDs.skullstripped
        WM_img = WM_img(startsWith(WM_img_name,'rc2s'));
    end
    matlabbatch{2}.spm.tools.dartel.warp.images{2}(i,1) = WM_img;
end

matlabbatch{2}.spm.tools.dartel.warp.settings.template = [IDs.templatename '_'  datestr(now,'yymmdd')];
matlabbatch{2}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{2}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{2}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{2}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{2}.spm.tools.dartel.warp.settings.optim.its = 3;

% purge variables
clear GM_img WM_img

%% SAVE AND RUN JOB
%------------------
save(fullfile(scriptpath,['DARTEL_' session{1} '_' session{3} '_' datestr(now,'yymmdd') '.mat']),'matlabbatch');
spm_jobman('run', matlabbatch, '' );
clear matlabbatch

end
