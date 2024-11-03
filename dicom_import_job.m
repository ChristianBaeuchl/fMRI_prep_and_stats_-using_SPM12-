% matlabbatch = dicom_import_job(raw,outdir)
% The function takes the list of raw files and output directory and passes
% them on to the matlabbatch sturcture array
%
% 

function matlabbatch = dicom_import_job(raw,outdir)

for i = 1:numel(outdir)    
    matlabbatch{1}.spm.util.dicom.data = raw;
    matlabbatch{1}.spm.util.dicom.root = 'flat';
    matlabbatch{1}.spm.util.dicom.outdir = outdir;
    matlabbatch{1}.spm.util.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;
end


