%% check images
% helps to automatize the consecutive visual inspection of co-registered
% and/or spatially normalized functional images
%
%  Christian Baeuchl | last edit 23.05.2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Closes all (plot) windows, clears variables and command window, starts to clock the scripts' running time 
close all
clear,clc
tStart = tic;

%% Define folder name, prefixes of functional and structural images and indicate
%  what kind of image check should be done:
%  - 'norm' : compare normalized images to each other
%  - 'coreg': check the alignment of the mean EPI to the structural image

FunctFolder = 'task folder'; 
FunctPrefix = '^swauf.*\.nii$'; %'mask.nii'; %'^uf.*\.nii$'
FunctMean   = '^mean.*\.nii$';
StructFolder = 'T1';
StructPrefix = '^s.*\.nii$';
ToCheck     = 'coreg'; %'norm'; %'coreg';

%% Path and IDs
datapath =  'path to your preprocessed data'; 
%IDs = cellstr(num2str(load('IDs_FOO_selected.txt')));
%IDs = cellstr(num2str((sort(importdata('IDs_test.txt')))));
IDs = (importdata('IDs_test.txt'));

if strcmp(ToCheck,'norm') || strcmp(ToCheck,'mask')
    ImgOnScreen = 8;
    NumOfScreens = floor(numel(IDs)/ImgOnScreen);
    remainder = rem(numel(IDs),ImgOnScreen);
end

%% Initialise SPM defaults
spm('defaults','fmri');
spm_jobman('initcfg');

%%
if strcmp(ToCheck,'coreg')
    
    for i = 1:numel(IDs)
        struct          = cellstr(spm_select('FPList',[datapath filesep IDs{i} '_1' filesep StructFolder], StructPrefix));
        funct_mean      = spm_select('FPList', fullfile(datapath,[IDs{i} '_1'],FunctFolder), FunctMean);
        funct_full_s2   = spm_select('FPList', fullfile(datapath,[IDs{i} '_2'],FunctFolder), FunctPrefix);
        funct_single_s2 = funct_full_s2(1,:);
        data            = cellstr([struct;funct_mean;funct_single_s2]);
        
        %%
        matlabbatch{1}.spm.util.checkreg.data = data;
        batch_matrix(i,:)                     = matlabbatch;
        clear struct funct_mean funct_full_s2 funct_single_s2 data matlabbatch
    end
    
    %%
    tEnd = toc(tStart);
    fprintf('\n ... done!\n');
    fprintf('\nProcess finished in %d hours, %d minutes, %d seconds and %d milliseconds\n',...
        floor(tEnd/(60*60)),floor(rem((tEnd/60),60)),floor(rem(tEnd,60)),round((tEnd-floor(tEnd))*1000));
    
    %%
    for j = 1:size(batch_matrix,1)
        
        fprintf('\nPlease press "enter" if you want to check the structural-to-functional co-registration of %s\n',IDs{j})
        pause;
        spm_jobman('run', batch_matrix(j,:), '' );
        fprintf('\n Window #%d/%d\n',j,size(batch_matrix,1))
        
    end
    
elseif strcmp(ToCheck,'norm') || strcmp(ToCheck,'mask')    
    n = 1;   
    if remainder == 0
        for i = 1:NumOfScreens
            for ii = 1:ImgOnScreen
                funct_fulllist   =  spm_select('FPList', fullfile(datapath,IDs{n},FunctFolder), FunctPrefix);
                funct_data(ii,1) = cellstr(funct_fulllist(1,:));
                n = n+1;
            end
            clear funct_fulllist
            
            %%
            matlabbatch{1}.spm.util.checkreg.data = funct_data;
            batch_matrix(i,:) = matlabbatch;
            clear funct_data matlabbatch
        end
        
    elseif remainder ~= 0
        for i = 1:(NumOfScreens+1)
            if i ~= (NumOfScreens+1)
                for ii = 1:ImgOnScreen
                    funct_fulllist   =  spm_select('FPList', fullfile(datapath,IDs{n},FunctFolder), FunctPrefix);
                    funct_data(ii,1) = cellstr(funct_fulllist(1,:));
                    n = n+1;
                end
            else
                for iii = 1:remainder
                    funct_fulllist   =  spm_select('FPList', fullfile(datapath,IDs{n},FunctFolder), FunctPrefix);
                    funct_data(iii,1) = cellstr(funct_fulllist(1,:));
                    n = n+1;
                end
            end
            clear funct_fulllist
            
            %%
            matlabbatch{1}.spm.util.checkreg.data = funct_data;
            batch_matrix(i,:) = matlabbatch;
            clear funct_data matlabbatch
        end
    end
    
    %%
    tEnd = toc(tStart);
    fprintf('\n ... done!\n');
    fprintf('\nProcess finished in %d hours, %d minutes, %d seconds and %d milliseconds\n',...
        floor(tEnd/(60*60)),floor(rem((tEnd/60),60)),floor(rem(tEnd,60)),round((tEnd-floor(tEnd))*1000));
    
    %%
    for j = 1:size(batch_matrix,1)
        
        if j == 1
            fprintf('\nPlease press "enter" if you want to check the first set of %d images\n',ImgOnScreen)
            pause;
            spm_jobman('run', batch_matrix(j,:), '' );
            fprintf('\n Window #%d\n',j)
        elseif j > 1 && j < size(batch_matrix,1)
            fprintf('Please press "enter" if you want to check the next set of %d images\n',ImgOnScreen)
            pause;
            spm_jobman('run', batch_matrix(j,:), '' );
            fprintf('\n Window #%d\n',j)
        else
            if remainder == 0
                fprintf('Please press "enter" if you want to check the last set of %d images\n',ImgOnScreen)
            elseif remainder ~= 0
                fprintf('Please press "enter" if you want to check the last set of %d images\n',remainder)
            end
            pause;
            spm_jobman('run', batch_matrix(j,:), '' );
            fprintf('\n Window #%d\n',j)
        end
        
    end
    
end

