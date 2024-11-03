%% Nifti 3D to 4D conversion script

% Christian Baeuchl - 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear variables and command window
clear; clc;
tStart = tic;

%% Specify paths
datapath   = 'path to your "Dicom imported" data';
scriptpath = 'path to your dicom_import_job script';
addpath(scriptpath)

X  = 3; % which machine? - For Windows: X = 3; For Unix: X = 1
TR = 2; % TR = time of repetition

folderCount = [];
Task  = 'task name';
Sess  = 'session name';
FPnum = 'number of funding periods to assess';


%%  Initialise SPM defaults
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%% loop over FPs
for i = 1:FPnum

    folderList = cellstr(num2str(sort(importdata(['.\Participant_lists\List_' Task '_' Sess '.txt']))));
    % folderList = {
    %     'subjcect_01'
    %     'subjcect_02'
    %     'subjcect_03'
    %     ...
    %     };

    for j = 1:length(folderList)

        scanList  = ls([datapath '0' num2str(i) '_*' '\01_Prep_data\' char(folderList{j})]);
        scanList  = scanList(X:end,:);
        toConvert = regexp(cellstr(scanList),'.*EPI$||.*State$','match');
        toConvert(cellfun(@isempty,toConvert)) = [];
        toConvert = toConvert(1); % ERASE LATER

        for k = 1:length(toConvert)

            threeDpath = fullfile(datapath,ls([datapath '0' num2str(i) '_*']),'01_Prep_data',char(folderList{j}),char(toConvert{k}),'3Dnifti');
            threeDdata = cellstr(spm_select('FPList',threeDpath,'.*nii$'));
            if isempty(char(threeDdata)) % jump to next iteration of loop if there are no data to convert for the current sequence
                continue
            end

            % display which subject and sequence is being processed
            if ~isempty(threeDdata)
                fprintf('Converting 3D niftis of %s data from subject "%s" in MR scanning session #%s (%s files)\n',...
                    char(toConvert{k}),char(folderList{j}),num2str(i),sprintf('%d',size(threeDdata,1)));
            end

            fourDpath      = fullfile(datapath,ls([datapath '0' num2str(i) '_*']),'01_Prep_data',char(folderList{j}),char(toConvert{k}));
            fourDname      = char(toConvert{k});
            indx           = strfind(fourDname,'_');
            fourDname      = fourDname(1:indx);
            pathAndName    = [fourDpath '\' char(folderList{j}) '_' fourDname '4D.nii'];
            matlabbatch(k) = threeDto4D_job(threeDdata,pathAndName,TR);

        end

        %% save and run job
        if any(~cellfun(@isempty,matlabbatch)) % remove empty cells
            matlabbatch = matlabbatch(~cellfun('isempty',matlabbatch));
        end
        save(fullfile(datapath,ls([datapath '0' num2str(i) '_*']),'01_Prep_data',char(folderList{j}),[char(folderList{j}) '_' extractAfter(ls([datapath '0' num2str(i) '_*']),'_') '_3Dto4D.mat']),'matlabbatch');
        spm_jobman('serial', matlabbatch);

        clear matlabbatch threeDpath threeDdata fourDpath fourDname pathAndName toConvert scanList

    end

    folderCount = [folderCount folderList];
    clear fulldirectory folderList

end

%%
cd(datapath)
tEnd = toc(tStart);
fprintf('\nFinished converting nifti images of %d subjects within %d hours, %d minutes, %d seconds and %d milliseconds\n',...
    size(folderCount,1),floor(tEnd/(60*60)),floor(rem((tEnd/60),60)),floor(rem(tEnd,60)),round((tEnd-floor(tEnd))*1000));
diary off

%% END
