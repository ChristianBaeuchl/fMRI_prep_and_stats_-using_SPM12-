%% findHolesIn3DconImgs
% this script takes binary images of individual 1st level SPM analyes
% and checks whether they have holes in them, using
% the Baycrest NIFTI toolbox and MATLAB's build-in function 'imfill'
%
% by Christian Baeuchl - 24.03.2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear , clc
tStart = tic;

%% Full path of experiment folder
path =  'path to your 1st level stats';
IDs = sort(importdata('IDs_list.txt'));

%% specify stats folder and image
statsFolder = '\1st level stats folder\';
Img         = 'mask.nii';

%% loop over IDs
for i = 1:numel(IDs)
    
    tmp00 = load_untouch_nii([path IDs{i} statsFolder Img]);
    tmp01 = tmp00;
    tmp01.img = imfill(tmp01.img);
    imgDiff = tmp01.img - tmp00.img;
    checkSum = sum(sum(sum(imgDiff)));
    
    if checkSum ~= 0
        fprintf('Subject/Session %s has at least one hole in his or her brain!\n',IDs{i})
    end
    
    % clear variables
    clear tmp00 tmp01 imgDiff checkSum
    
end

tEnd = toc(tStart);
fprintf('\n ... done!\n');