%% DICOM Import script
%  This script converts DICOM files into nifti Files using SPMs DICOM
%  Import sequence. It follows a particular directory structure which must
%  be adhered to.

% Christian Baeuchl - 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear variables and command window
clear; clc;
tStart = tic;

%% Specify paths
% In- and outputfolder for RAW data
FP = 2; % 2 , 3. Which funding period?
datapath   = 'path to your dicom images'; 
NIC_ID     = 'folder with project ID';
outpath    = ['path that data are exported to _' num2str(FP) '- uses funding period information in path name'];
scriptpath = 'path of dicom_import_job script'; 
addpath(scriptpath)

% Shall data be imported for subject who already have a folder in outpath?
reimport = 1; % 0 = don't import data again; % 1 = import data again

% Create a diary of command window outputs and store it in the output path
diary([outpath filesep 'dicom_import_diary_' datestr(now,'yymmdd') '.txt']);

%% Which OS do you run this script on? (relevant for 'ls' function)
X = 3; % For Windows: X = 3; For Unix: X = 1

%%  Subject list (and respective folder)
% Either 
%IDs = cellstr(num2str(load('List_with_subject_folders.txt')));
IDs = {
    'subject_01'
    'subject_02'
    'subject_03'
    ...
    };

%%  Sequence folders
sequences = {
    't1_mpr_ns_sag_iso_m'                        % MP-RAGE
    's0[1-9]\w*gre_field_3x3x3_2g0_8_38sl_m'     % fieldmap magnitude image
    's0[1-9]\w*gre_field_3x3x3_2g0_8_38sl_p'     % fieldmap phase image
    'rs_epi_3x3x3_2g0_8_34sl_m'                  % resting-state EPI
    'gng_epi_3x3x3_2g0_8_34sl_m'                 % go/nogo task EPI
    'str_isi_epi_3x3x3_2g0_8_34sl_m'             % stroop task EPI
    'sst_isi_epi_3x3x3_2g0_8_34sl_m'             % stop-signal task EPI
    'hare_epi_3x3x3_2g0_8_34sl_m'                % complex decision making task EPI
    's[1-2][1-9]\w*gre_field_3x3x3_2g0_8_38sl_m' % fieldmap magnitude image #2
    's[1-2][1-9]\w*gre_field_3x3x3_2g0_8_38sl_p' % fieldmap phase image #2
    's15_mbep2d_diff_AP_1_7iso_67dir_b1000_d'    % the rest are DTI images
    's16_mbep2d_diff_AP_1_7iso_67dir_b1000_d'
    's17_mbep2d_diff_AP_1_7iso_67dir_b1000_d'
    's18_mbep2d_diff_AP_1_7iso_67dir_b1000_d'
    's19_mbep2d_diff_AP_1_7iso_67dir_b1000_d'
    's20_mbep2d_diff_AP_1_7iso_67dir_b1000_d'
    'mbep2d_diff_PA_1_7iso_6dir_b0_d'
    };

dti_idx = 11:16;

sequences_folder = {
    'T1'
    'Fieldmap\Magnitude'
    'Fieldmap\Phase'
    'Resting_State'
    'GoNogo_EPI'
    'Stroop_EPI'
    'StopSignal_EPI'
    'ComplexDM_EPI'
    'Fieldmap_2nd\Magnitude'
    'Fieldmap_2nd\Phase'
    'Diff_Weighted\AP_scan_01'
    'Diff_Weighted\AP_scan_02'
    'Diff_Weighted\AP_scan_03'
    'Diff_Weighted\AP_scan_04'
    'Diff_Weighted\AP_scan_05'
    'Diff_Weighted\AP_scan_06'
    'Diff_Weighted\PA_scan'
    };

del_imgs = 0; % 0 if you do not want to delete images after import

%%  Initialise SPM defaults
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% initialize matlabbatch cells
matlabbatch  = cell(1,size(sequences,1));
batch_matrix = cell(numel(IDs),size(sequences,1));

%% loop over IDs
for i = 1:length(IDs)

    % skip current subject if folder already exists in outpath
    if reimport == 0
        type = exist(fullfile(outpath,IDs{i}),'dir');
        if type == 7
            continue
        end
    end

    D = dir([datapath filesep NIC_ID '*' filesep '*' IDs{i} '*']);


    if size(D,1) == 1 && FP == 1
        k = 1;
    elseif size(D,1) > 1 && FP == 1
        k = find(~cellfun(@isempty,regexp({D.name},'^211_*\d{5}$','match','once')));
        if isempty(k)
            k = find(~cellfun(@isempty,regexp({D.name},'^211_*\d{5}_1$','match','once')));
        end
        if length(k) == 2 && (datenum(D(1).date) < datenum(D(2).date))
            k = 1;
        elseif length(k) == 2 && (datenum(D(1).date) > datenum(D(2).date))
            k = 2;
        elseif length(k) > 2
            fprintf(2,'Check folder %s manually - skipping data set!\n',D(1).name)
            k = [];
            continue
        end
    elseif size(D,1) == 1 && FP == 2
        k = 1;
    elseif size(D,1) > 1 && FP == 2
        k = find(~cellfun(@isempty,regexp({D.name},'^211_*\d{5}_2$','match','once')));
        if isempty(k)
            k = find(~cellfun(@isempty,regexp({D.name},'^211_*\d{5}$','match','once')));
        end
        if length(k) == 2 && (datenum(D(1).date) > datenum(D(2).date))
            k = 1;
        elseif length(k) == 2 && (datenum(D(1).date) < datenum(D(2).date))
            k = 2;
        elseif k > 2
            fprintf(2,'Check folder %s manually - skipping data set!\n',D(1).name)
            k = [];
            continue
        end
    elseif size(D,1) > 1 && FP == 3
        k = find(~cellfun(@isempty,regexp({D.name},'^211_*\d{5}_3$','match','once')));
        if isempty(k)
            k = find(~cellfun(@isempty,regexp({D.name},'^211_*\d{5}_[4-9]$','match','once')));
            if ~(numel(k) == 1 && year(datetime(D(k).date)) >= 2020)
                fprintf(2,'Check folder %s manually - skipping data set!\n',D(1).name)
                k = [];
                continue
            end
        end
        if isempty(k)
            k = find(~cellfun(@isempty,regexp({D.name},'^211_*\d{5}$','match','once')));
        end
    end

    if ~exist('k','var')
        fprintf(2,'Warning - folder name does not match subject ID in funding period %d!\n',FP)
        fprintf(2,'Skipping data set of ID %s!\n',IDs{i})
        continue
    end

    % skip unnecessary folders
    wd = cd([D(k).folder filesep D(k).name]);
    datefolder = ls;
    datefolder = datefolder(X:end,:);

    if size(datefolder,1) == 1 && FP == 1
        if str2num(datefolder(1,1:4)) > 2016
            fprintf(2,'Warning - there is only one data set for %s, but it is not from FP1 %d according to its folder date!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            datefolder_selected = datefolder;
        end
    elseif size(datefolder,1) > 1 && size(datefolder,1) == 2 && FP == 1
        if strcmp(datefolder(1,1:4),datefolder(2,1:4)) % in case both date folders should start with the same year
            fprintf(2,'Warning - there are 2 date folders for %s for FP %d that indicate the same year!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            datefolder_selected = num2str(min(str2num(datefolder)));
        end
    elseif size(datefolder,1) > 1 && size(datefolder,1) == 3 && FP == 1
        if strcmp(datefolder(1,1:4),datefolder(2,1:4)) ||  strcmp(datefolder(1,1:4),datefolder(3,1:4)) ...
                || strcmp(datefolder(2,1:4),datefolder(3,1:4)) % in case any of the 3 date folders should start with the same year
            fprintf(2,'Warning - there are 3 date folders for %s for FP %d of which at least 2 indicate the same year!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            datefolder_selected = num2str(min(str2num(datefolder)));
        end
    elseif size(datefolder,1) == 1 && FP == 2
        if str2num(datefolder(1,1:4)) < 2016 || str2num(datefolder(1,1:4)) > 2020
            fprintf(2,'Warning - there is only one data set for %s, but it is not from FP %d according to its folder date!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            datefolder_selected = datefolder;
        end
    elseif size(datefolder,1) > 1 && size(datefolder,1) == 2 && FP == 2
        if strcmp(datefolder(1,1:4),datefolder(2,1:4)) % in case both date folders should start with the same year
            fprintf(2,'Warning - there are 2 date folders for %s for FP %d that indicate the same year!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            folderyear1 =  (str2num(datefolder(1,1:4)) >= 2016 && str2num(datefolder(1,1:4)) <= 2020);
            folderyear2 =  (str2num(datefolder(2,1:4)) >= 2016 && str2num(datefolder(2,1:4)) <= 2020);
            if folderyear1 && folderyear2
                fprintf(2,'Warning, there are 2 folders who both fall into the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            elseif folderyear1
                datefolder_selected = datefolder(1,:);
            elseif folderyear2
                datefolder_selected = datefolder(2,:);
            else
                fprintf(2,'Warning both folders are not within the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            end
        end
    elseif size(datefolder,1) > 1 && size(datefolder,1) == 3 && FP == 2
        if strcmp(datefolder(1,1:4),datefolder(2,1:4)) ||  strcmp(datefolder(1,1:4),datefolder(3,1:4)) ...
                || strcmp(datefolder(2,1:4),datefolder(3,1:4)) % in case any of the 3 date folders should start with the same year
            fprintf(2,'Warning - there are 3 date folders for %s for FP %d of which at least 2 indicate the same year!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            folderyear1 =  (str2num(datefolder(1,1:4)) >= 2016 && str2num(datefolder(1,1:4)) <= 2020);
            folderyear2 =  (str2num(datefolder(2,1:4)) >= 2016 && str2num(datefolder(2,1:4)) <= 2020);
            folderyear3 =  (str2num(datefolder(3,1:4)) >= 2016 && str2num(datefolder(3,1:4)) <= 2020);
            if (folderyear1 && folderyear2 && folderyear3) || (folderyear1 && folderyear2) || (folderyear1 && folderyear3)...
                    || (folderyear2 && folderyear3)
                fprintf(2,'Warning, there are at least 2 folders who both fall into the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            elseif folderyear1
                datefolder_selected = datefolder(1,:);
            elseif folderyear2
                datefolder_selected = datefolder(2,:);
            elseif folderyear3
                datefolder_selected = datefolder(3,:);
            else
                fprintf(2,'Warning one of the folder are within the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            end
        end
    elseif size(datefolder,1) == 1 && FP == 3
        if str2num(datefolder(1,1:4)) < 2020
            fprintf(2,'Warning - there is only one data set for %s, but it is not from FP %d according to its folder date!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            datefolder_selected = datefolder;
        end
    elseif size(datefolder,1) > 1 && size(datefolder,1) == 2 && FP == 3
        if strcmp(datefolder(1,1:4),datefolder(2,1:4)) % in case both date folders should start with the same year
            fprintf(2,'Warning - there are 2 date folders for %s for FP %d that indicate the same year!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            folderyear1 =  (str2num(datefolder(1,1:4)) >= 2020 && str2num(datefolder(1,1:4)) <= 2024);
            folderyear2 =  (str2num(datefolder(2,1:4)) >= 2020 && str2num(datefolder(2,1:4)) <= 2024);
            if folderyear1 && folderyear2
                fprintf(2,'Warning, there are 2 folders who both fall into the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            elseif folderyear1
                datefolder_selected = datefolder(1,:);
            elseif folderyear2
                datefolder_selected = datefolder(2,:);
            else
                fprintf(2,'Warning both folders are not within the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            end
        end
    elseif size(datefolder,1) > 1 && size(datefolder,1) == 3 && FP == 3
        if strcmp(datefolder(1,1:4),datefolder(2,1:4)) ||  strcmp(datefolder(1,1:4),datefolder(3,1:4)) ...
                || strcmp(datefolder(2,1:4),datefolder(3,1:4)) % in case any of the 3 date folders should start with the same year
            fprintf(2,'Warning - there are 3 date folders for %s for FP %d of which at least 2 indicate the same year!\n',IDs{i},FP);
            fprintf(2,'Skipping data set\n');
            k = [];
            continue
        else
            folderyear1 =  (str2num(datefolder(1,1:4)) >= 2020 && str2num(datefolder(1,1:4)) <= 2024);
            folderyear2 =  (str2num(datefolder(2,1:4)) >= 2020 && str2num(datefolder(2,1:4)) <= 2024);
            folderyear3 =  (str2num(datefolder(3,1:4)) >= 2020 && str2num(datefolder(3,1:4)) <= 2024);
            if (folderyear1 && folderyear2 && folderyear3) || (folderyear1 && folderyear2) || (folderyear1 && folderyear3)...
                    || (folderyear2 && folderyear3)
                fprintf(2,'Warning, there are at least 2 folders who both fall into the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            elseif folderyear1
                datefolder_selected = datefolder(1,:);
            elseif folderyear2
                datefolder_selected = datefolder(2,:);
            elseif folderyear3
                datefolder_selected = datefolder(3,:);
            elseif size(datefolder,1) > 3
                fprintf(2,'Warning - there are more than 3 date folders for %s !\n',IDs{i});
                fprintf(2,'Please check manually\n');
                fprintf(2,'Skipping data set\n')
                k = [];
                continue
            else
                fprintf(2,'Warning all 3 folders are not within the year range of FP %d!\n',FP);
                fprintf(2,'Skipping data set\n');
                k = [];
                continue
            end
        end
    end


    wd = cd([D(k).folder filesep D(k).name filesep datefolder_selected]);
    dummyfolder = ls;

    %% iterate over several subfolders if data files are split up
    for x = X:size(dummyfolder,1)

        wd = cd([D(k).folder filesep D(k).name filesep datefolder_selected  filesep dummyfolder(x,:)]);
        wd = pwd;
        disp(wd);

        % get list of folders of current subject & current session
        LS_ID = ls;

        %% load raw data for each sequence without a prompt
        for j = 1:numel(sequences)

            % change into folder of current subject & current session
            cd(wd)

            n = 1;
            for h = X:size(LS_ID,1)
                if ismember(j,dti_idx)
                    break
                else
                    % test = strfind(LS_ID(h,:),sequences{j});
                    test = regexp(LS_ID(h,:),sequences{j},'match');
                    if ~isempty(test)
                        curr_seq(n,:) =  LS_ID(h,:); % grabs the current sequence file
                        n = n+1;
                    end
                end
            end

            z = 0;
            if FP == 2 || FP == 3
                for g = X:size(LS_ID,1)
                    if ~ismember(j,dti_idx)
                        break
                    else
                        if z == 0
                            test2 = regexp(LS_ID(g,:),sequences{j},'match');
                            if ~isempty(test2)
                                dti_num = str2num(LS_ID(g,2:3))+1;
                                curr_seq(n,:) = LS_ID(g,:); % grabs the current sequence file
                                n = n+1;
                                z = z+1;
                                test2 = [];
                            end
                        elseif z ~= 0
                            test2 = regexp(LS_ID(g,:),['s' num2str(dti_num) '_' sequences{j}] ,'match');
                            if ~isempty(test2)
                                curr_seq(n,:) = LS_ID(g,:); % grabs the current sequence file
                                dti_num = dti_num+1;
                                n = n+1;
                            end
                        end
                    end
                end
            end

            if ~exist('curr_seq','var') % jumps to next iteration of the loop if one of the MRI sequences cannot be found in this directory
                continue
            end

            % skip unnecessary folders
            if size(curr_seq,1) == 1
                wd_new = cd([wd filesep curr_seq]);
                wd_new = pwd;
                disp(wd_new);
                dummyfolder2 = ls;
                % get the raw data
                raw = cellstr(spm_select('FPList',fullfile(wd_new,dummyfolder2(X,:)),'0*\.ima*')); % lists all files in indicated path that start with a '0' and end with .ima
                check_dcm = dicominfo(raw{1});
                if ~strcmp(check_dcm.PatientName.GivenName,IDs{i})
                    fprintf(2,'ID of import data (%s) does not match subject ID (%s)!\n',check_dcm.PatientName.GivenName,IDs{i});
                    fprintf(2,'Skipping dataset\n');
                    k = [];
                    continue
                end
            elseif size(curr_seq,1) > 1
                warning('Multiple folders found for the same sequence:')
                fprintf('\n')
                for m = 1:size(curr_seq,1)
                    fprintf('%s\n',curr_seq(m,:))
                end
                fprintf('\nTaking data from folder with largest amount of files\n')
                fprintf(2,'\nPlease check if this is actually the folder whose content is supposed to be imported!!!\n')
                for m = 1:size(curr_seq,1)
                    wd_new = cd([wd filesep curr_seq(m,:)]);
                    wd_new = pwd;
                    disp(wd_new);
                    dummyfolder4 = ls;
                    sizefolder(m) = size(cellstr(spm_select('FPList',fullfile(wd_new,dummyfolder4(X,:)),'0*\.ima*')),1);
                end
                wd_new = fileparts(wd_new);
                if sum(sizefolder == max(sizefolder)) > 1
                    ind = find(sizefolder == max(sizefolder));
                    wd_new = cd([wd_new filesep curr_seq(ind(1),:)]);
                elseif sum(sizefolder == max(sizefolder)) == 1
                    wd_new = cd([wd_new filesep curr_seq(sizefolder == max(sizefolder),:)]);
                end
                wd_new = pwd;
                disp(wd_new);
                dummyfolder5 = ls;
                raw = cellstr(spm_select('FPList',fullfile(wd_new,dummyfolder5(X,:)),'0*\.ima*'));
                check_dcm = dicominfo(raw{1});
                if ~strcmp(check_dcm.PatientName.GivenName,IDs{i})
                    fprintf(2,'ID of import data (%s) does not match subject ID (%s)!\n',check_dcm.PatientName.GivenName,IDs{i});
                    fprintf(2,'Skipping dataset\n');
                    k = [];
                    continue
                end
            end

            if isempty(raw)
                fprintf(2,'Warning - Could not grab DICOM files for subject $s, FP %d, sequence %s\n',IDs{i},FP,sequence{j})
            end

            % select appropriate output folder depending on the MRI sequence
            if (~strcmp(sequences_folder(j),'T1') && ~strcmp(sequences_folder(j),'Fieldmap\Magnitude')...
                    && ~strcmp(sequences_folder(j),'Fieldmap\Phase') && ~strcmp(sequences_folder(j),'Fieldmap_2nd\Magnitude')...
                    && ~strcmp(sequences_folder(j),'Fieldmap_2nd\Phase') && ~strcmp(sequences_folder(j),'Diff_Weighted\AP_scan_01')...
                    && ~strcmp(sequences_folder(j),'Diff_Weighted\AP_scan_01') && ~strcmp(sequences_folder(j),'Diff_Weighted\AP_scan_02')...
                    && ~strcmp(sequences_folder(j),'Diff_Weighted\AP_scan_03') && ~strcmp(sequences_folder(j),'Diff_Weighted\AP_scan_04')...
                    && ~strcmp(sequences_folder(j),'Diff_Weighted\AP_scan_05') && ~strcmp(sequences_folder(j),'Diff_Weighted\AP_scan_06')...
                    && ~strcmp(sequences_folder(j),'Diff_Weighted\PA_scan'))
                type = exist(fullfile(outpath,IDs{i},char(sequences_folder(j)),'3Dnifti'),'dir');
                if type ~= 7
                    mkdir(fullfile(outpath,IDs{i},char(sequences_folder(j)),'3Dnifti'));
                end
                outdir = cellstr(fullfile(outpath,IDs{i},sequences_folder(j),'3Dnifti'));
            else
                type = exist(fullfile(outpath,IDs{i},char(sequences_folder(j))),'dir');
                if type ~= 7
                    mkdir(fullfile(outpath,IDs{i},char(sequences_folder(j))));
                end
                outdir = cellstr(fullfile(outpath,IDs{i},sequences_folder(j)));
            end


            % display which subject and sequence is being processed
            if size(curr_seq,1) == 1
                fprintf('Importing subject "%s", "%s" (%s files)\n',...
                    IDs{i},curr_seq,sprintf('%d',size(raw,1)));

            elseif size(curr_seq,1) > 1
                fprintf('Importing subject "%s", "%s" (%s files)\n',...
                    IDs{i},curr_seq(sizefolder == max(sizefolder),:),sprintf('%d',size(raw,1)));
            end

            % clear variables to avoid errors
            clear curr_seq sizefolder

            %% CREATE MATLABBATCH
            matlabbatch(:,j) = dicom_import_job(raw,outdir);
        end
        % clear variables to avoid errors
        clear raw
    end
    % clear variables to avoid errors
    clear D k

    %% SAVE JOB
    %----------------------------------------------------------------------
    if size(matlabbatch,2) ~= size(sequences,1)
        matlabbatch(1,size(matlabbatch,2)+1:size(sequences,1)) = cell(1,1);
    end
    batch_matrix(i,:) = matlabbatch;
    if all(~cellfun(@isempty,matlabbatch))
        save(fullfile(outpath,IDs{i},'dicom_import.mat'),'matlabbatch');
    end
    matlabbatch = cell(1,1);

end

%% Create list of missing data sets
empty_rows  = find(all(cellfun(@isempty,batch_matrix),2));
[rows,cols] = find(cellfun(@isempty,batch_matrix));
[idx,test]  = ismember(rows,empty_rows);
part_empty_rows      = rows;
part_empty_rows(idx) = [];

part_uniq = unique(part_empty_rows);
part_idx = zeros(size(IDs,1),1);
for k=1:numel(part_uniq)
    part_m(k) = numel(part_empty_rows(part_empty_rows==part_uniq(k)));
    part_idx(part_uniq(k)) = 1;
end

missing_data   = IDs(all(cellfun(@isempty,batch_matrix),2));
missing_matrix = zeros(size(missing_data,1),size(sequences,1));
missing_part   = IDs(logical(part_idx));
missing_data   = [missing_data;missing_part];
part_matrix    = zeros(size(missing_part,1),size(sequences,1));
for k=1:size(part_matrix,1)
    part_matrix(k,1:(numel(sequences)-part_m(k))) = 1;
end
missing_matrix = [missing_matrix;part_matrix];

tableMM = array2table(missing_matrix);
tableMM.Properties.VariableNames = sequences;
tableID = table(missing_data);
tableID.Properties.VariableNames = {'IDs'};
tableMM = [tableID,tableMM];
tableMM = sortrows(tableMM);
writetable(tableMM,[outpath filesep 'FP_' num2str(FP) '_missing_data_list_' datestr(now,'yymmdd_HHMM') '.xlsx']);
disp(tableMM)

%% Delete rows will all empty cells
batch_matrix(all(cellfun(@isempty,batch_matrix),2), : ) = [];


%% RUN JOB
%------------------

parfor i = 1:size(batch_matrix,1)
    try
        spm_jobman('serial', batch_matrix(i,:), '' );
    catch
        sprintf(['Could not run DICOM import for subject: ' IDs{i}]);
    end
end

if del_imgs
    for i = 1:size(batch_matrix,1)
        for d = 1:del_imgs
            delete(tmp(d,:));
        end
    end
end

%%
cd(outpath)
tEnd = toc(tStart);
fprintf('\nFinished importing DICOM images of %d subjects within %d hours, %d minutes, %d seconds and %d milliseconds\n',...
    size(IDs,1),floor(tEnd/(60*60)),floor(rem((tEnd/60),60)),floor(rem(tEnd,60)),round((tEnd-floor(tEnd))*1000));
diary off
%% END
