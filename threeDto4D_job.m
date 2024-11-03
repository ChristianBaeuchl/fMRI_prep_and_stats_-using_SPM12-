%-----------------------------------------------------------------------
% Job saved on 23-Jun-2023 11:31:33 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

function matlabbatch = threeDto4D_job(vols,PathAndName,TR)

matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = PathAndName;
matlabbatch{1}.spm.util.cat.dtype = 0;
matlabbatch{1}.spm.util.cat.RT    = TR;

end