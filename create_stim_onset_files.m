%% Create stimulus onset files
%
% This script creates stimulus onset files for 4 fMRI tasks.
% It takes event onset information from the output of the 4 "analyze"
% scripts, namely analyze_CDM, analyze_GoNogo, analyze_StopSignal and
% analyze_Stroop. Be aware that errors in those outputs will naturally
% transfer to errors in the files created by this script.
%
% Christian Baeuchl - 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear variables, close windows and clear the command line
clear,close, clc

%% Define paths and load data
datapath = 'path to the code that preprocessed your taks';
outpath  = 'path to where your stimulus onset files should end up';

CDM_BL  = load([datapath 'CDM_Onsets_BL.mat']);
CDM_FU3 = load([datapath 'CDM_Onsets_FU3.mat']);
CDM_FU6 = load([datapath 'CDM_Onsets_FU6.mat']);

GoNogo_BL  = load([datapath 'GoNogo_Onsets_BL.mat']);
GoNogo_FU3 = load([datapath 'GoNogo_Onsets_FU3.mat']);
GoNogo_FU6 = load([datapath 'GoNogo_Onsets_FU6.mat']);

Stroop_BL  = load([datapath 'Stroop_Onsets_BL.mat']);
Stroop_FU3 = load([datapath 'Stroop_Onsets_FU3.mat']);
Stroop_FU6 = load([datapath 'Stroop_Onsets_FU6.mat']);

StopSignal_BL  = load([datapath 'StopSignal_Onsets_BL.mat']);
StopSignal_FU3 = load([datapath 'StopSignal_Onsets_FU3.mat']);
StopSignal_FU6 = load([datapath 'StopSignal_Onsets_FU6.mat']);

% get rid of empty rows in the arrays of task data structures
fun = @(x) all(structfun(@isempty,x));

CDM_BL.CDM_Onsets(:,arrayfun(fun,CDM_BL.CDM_Onsets))   = [];
CDM_FU3.CDM_Onsets(:,arrayfun(fun,CDM_FU3.CDM_Onsets)) = [];
CDM_FU6.CDM_Onsets(:,arrayfun(fun,CDM_FU6.CDM_Onsets)) = [];

GoNogo_BL.GoNogo_Onsets(:,arrayfun(fun,GoNogo_BL.GoNogo_Onsets))   = [];
GoNogo_FU3.GoNogo_Onsets(:,arrayfun(fun,GoNogo_FU3.GoNogo_Onsets)) = [];
GoNogo_FU6.GoNogo_Onsets(:,arrayfun(fun,GoNogo_FU6.GoNogo_Onsets)) = [];

Stroop_BL.Stroop_Onsets(:,arrayfun(fun,Stroop_BL.Stroop_Onsets))   = [];
Stroop_FU3.Stroop_Onsets(:,arrayfun(fun,Stroop_FU3.Stroop_Onsets)) = [];
Stroop_FU6.Stroop_Onsets(:,arrayfun(fun,Stroop_FU6.Stroop_Onsets)) = [];

StopSignal_BL.StopSignal_Onsets(:,arrayfun(fun,StopSignal_BL.StopSignal_Onsets))   = [];
StopSignal_FU3.StopSignal_Onsets(:,arrayfun(fun,StopSignal_FU3.StopSignal_Onsets)) = [];
StopSignal_FU6.StopSignal_Onsets(:,arrayfun(fun,StopSignal_FU6.StopSignal_Onsets)) = [];

%
sessTasks = sort(who('CDM_BL','CDM_FU3','CDM_FU6','GoNogo_BL','GoNogo_FU3',...
    'GoNogo_FU6','Stroop_BL','Stroop_FU3','Stroop_FU6','StopSignal_BL',...
    'StopSignal_FU3','StopSignal_FU6'));

n = 1; % Initialize index for missing regressor messages
m = 1; % Initialize index for missing data set messages

%%
sepSess = any(contains(sessTasks,'_BL')) + any(contains(sessTasks,'_FU3')) + any(contains(sessTasks,'_FU6'));
numSess = numel(sessTasks);
fprintf(1,'\nProcessing stimulus onset files of %d sessions\n',sepSess)
fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  \n');
fprintf(1,'%d Task/Session combinations in total\n',numSess)
fprintf('\n')
fprintf('\n')

%% Loop over sessions
for i = 1:numSess

    if strcmp(sessTasks{i},'CDM_BL')

        IDs = [CDM_BL.CDM_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Complex Decision Making task (Baseline session)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %%
            if isempty(CDM_BL.CDM_Onsets(j).SubjectID)
                continue % skip iteration if data set for that is missing
            end

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2020 | 1st GLM

            %%
            names{1}         = 'Stimulus Decision';
            onsets{1}        = CDM_BL.CDM_Onsets(j).StimulusDecision;
            durations{1}     = CDM_BL.CDM_Onsets(j).RT_StimulusDecision;
            pmod(1).name{1}  = 'DecisionValue';
            pmod(1).param{1} = CDM_BL.CDM_Onsets(j).StimulusDecisionParametricReversed;
            pmod(1).poly{1}  = 1;
            orth{1}          = true;

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).MissedTrials))
                names{2}     = 'Missed Trials';
                onsets{2}    = CDM_BL.CDM_Onsets(j).MissedTrials;
                durations{2} = zeros(numel(CDM_BL.CDM_Onsets(j).MissedTrials),1);
                orth{2}      = true;
            end

            %% save onset files
            save([outpath num2str(CDM_BL.CDM_Onsets(j).SubjectID) '_CDM_BL_mk2020_GLM01' '.mat'], 'names','onsets', 'durations','pmod','orth');
            %% purge variables
            clear pmod

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2020 | 2nd GLM
            pmod(1).name{1}  = 'StimulusRatingReversed.ShortTerm';
            pmod(1).param{1} = CDM_BL.CDM_Onsets(j).StimulusDecisionParametricRatingsReversed(:,1);
            pmod(1).poly{1}  = 1;
            pmod(1).name{2}  = 'StimulusRatingReversed.LongTerm';
            pmod(1).param{2} = CDM_BL.CDM_Onsets(j).StimulusDecisionParametricRatingsReversed(:,2);
            pmod(1).poly{2}  = 1;

            %% save onset files
            save([outpath num2str(CDM_BL.CDM_Onsets(j).SubjectID) '_CDM_BL_mk2020_GLM02' '.mat'], 'names','onsets', 'durations','pmod','orth');
            %% purge variables
            clear names onsets durations pmod orth

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2021

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).SelfControlSuccessLongterm))
                names{1}     = 'Self Control Success | Longterm';
                onsets{1}    = CDM_BL.CDM_Onsets(j).SelfControlSuccessLongterm;
                durations{1} = zeros(numel(CDM_BL.CDM_Onsets(j).SelfControlSuccessLongterm),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: SelfControlSuccessLongterm | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).SelfControlSuccessShortterm))
                names{2}     = 'Self Control Success | Shortterm';
                onsets{2}    = CDM_BL.CDM_Onsets(j).SelfControlSuccessShortterm;
                durations{2} = zeros(numel(CDM_BL.CDM_Onsets(j).SelfControlSuccessShortterm),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: SelfControlSuccessShortterm | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).SelfControlFailureLongterm))
                names{3}     = 'Self Control Failure | Longterm';
                onsets{3}    = CDM_BL.CDM_Onsets(j).SelfControlFailureLongterm;
                durations{3} = zeros(numel(CDM_BL.CDM_Onsets(j).SelfControlFailureLongterm),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: SelfControlFailureLongterm | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).SelfControlFailureShortterm))
                names{4}     = 'Self Control Failure | Shortterm';
                onsets{4}    = CDM_BL.CDM_Onsets(j).SelfControlFailureShortterm;
                durations{4} = zeros(numel(CDM_BL.CDM_Onsets(j).SelfControlFailureShortterm),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: SelfControlFailureShortterm | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).NoConflictAllPositiveYes))
                names{5}     = 'No Conflict | All Positive | Accepted';
                onsets{5}    = CDM_BL.CDM_Onsets(j).NoConflictAllPositiveYes;
                durations{5} = zeros(numel(CDM_BL.CDM_Onsets(j).NoConflictAllPositiveYes),1);
                orth{5}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: NoConflictAllPositiveYes | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).NoConflictAllPositiveNo))
                names{6}     = 'No Conflict | All Positive | Declined';
                onsets{6}    = CDM_BL.CDM_Onsets(j).NoConflictAllPositiveNo;
                durations{6} = zeros(numel(CDM_BL.CDM_Onsets(j).NoConflictAllPositiveNo),1);
                orth{6}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: NoConflictAllPositiveNo | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).NoConflictAllNegativeYes))
                names{7}     = 'No Conflict | All Negative | Accepted';
                onsets{7}    = CDM_BL.CDM_Onsets(j).NoConflictAllNegativeYes;
                durations{7} = zeros(numel(CDM_BL.CDM_Onsets(j).NoConflictAllNegativeYes),1);
                orth{7}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: NoConflictAllNegativeYes | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_BL.CDM_Onsets(j).NoConflictAllNegativeNo))
                names{8}     = 'No Conflict | All Negative | Declined';
                onsets{8}    = CDM_BL.CDM_Onsets(j).NoConflictAllNegativeNo;
                durations{8} = zeros(numel(CDM_BL.CDM_Onsets(j).NoConflictAllNegativeNo),1);
                orth{8}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | BL | No onset of: NoConflictAllNegativeNo | Subject: %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            if exist('onsets','var')
                %% Remove empty cells (if they exist)
                names     = names(~cellfun('isempty',names));
                onsets    = onsets(~cellfun('isempty',onsets));
                durations = durations(~cellfun('isempty',durations));
                orth      = orth(~cellfun('isempty',orth));
                %% save onset files
                save([outpath num2str(CDM_BL.CDM_Onsets(j).SubjectID) '_CDM_BL_mk2021' '.mat'], 'names','onsets', 'durations','orth');
                %% purge variables
                clear names onsets durations orth
            else
                MissingDataSet{m,1} = sprintf('Due to missing data, cannot create onset file | _CDM_BL_mk2021 | %s', num2str(CDM_BL.CDM_Onsets(j).SubjectID));
                m = m+1;
            end

        end


    elseif strcmp(sessTasks{i},'CDM_FU3')

        IDs = [CDM_FU3.CDM_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Complex Decision Making task (Third follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %%
            if isempty(CDM_FU3.CDM_Onsets(j).SubjectID)
                continue % skip iteration if data set for that is missing
            end

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2020 | 1st GLM

            %%
            names{1}         = 'Stimulus Decision';
            onsets{1}        = CDM_FU3.CDM_Onsets(j).StimulusDecision;
            durations{1}     = CDM_FU3.CDM_Onsets(j).RT_StimulusDecision;
            pmod(1).name{1}  = 'DecisionValue';
            pmod(1).param{1} = CDM_FU3.CDM_Onsets(j).StimulusDecisionParametricReversed;
            pmod(1).poly{1}  = 1;
            orth{1}          = true;

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).MissedTrials))
                names{2}     = 'Missed Trials';
                onsets{2}    = CDM_FU3.CDM_Onsets(j).MissedTrials;
                durations{2} = zeros(numel(CDM_FU3.CDM_Onsets(j).MissedTrials),1);
                orth{2}      = true;
            end

            %% save onset files
            save([outpath num2str(CDM_FU3.CDM_Onsets(j).SubjectID) '_CDM_FU3_mk2020_GLM01' '.mat'], 'names','onsets', 'durations','pmod','orth');
            %% purge variables
            clear pmod

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2020 | 2nd GLM
            pmod(1).name{1}  = 'StimulusRatingReversed.ShortTerm';
            pmod(1).param{1} = CDM_FU3.CDM_Onsets(j).StimulusDecisionParametricRatingsReversed(:,1);
            pmod(1).poly{1}  = 1;
            pmod(1).name{2}  = 'StimulusRatingReversed.LongTerm';
            pmod(1).param{2} = CDM_FU3.CDM_Onsets(j).StimulusDecisionParametricRatingsReversed(:,2);
            pmod(1).poly{2}  = 1;

            %% save onset files
            save([outpath num2str(CDM_FU3.CDM_Onsets(j).SubjectID) '_CDM_FU3_mk2020_GLM02' '.mat'], 'names','onsets', 'durations','pmod','orth');
            %% purge variables
            clear names onsets durations pmod orth

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2021

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).SelfControlSuccessLongterm))
                names{1}     = 'Self Control Success | Longterm';
                onsets{1}    = CDM_FU3.CDM_Onsets(j).SelfControlSuccessLongterm;
                durations{1} = zeros(numel(CDM_FU3.CDM_Onsets(j).SelfControlSuccessLongterm),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: SelfControlSuccessLongterm | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).SelfControlSuccessShortterm))
                names{2}     = 'Self Control Success | Shortterm';
                onsets{2}    = CDM_FU3.CDM_Onsets(j).SelfControlSuccessShortterm;
                durations{2} = zeros(numel(CDM_FU3.CDM_Onsets(j).SelfControlSuccessShortterm),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: SelfControlSuccessShortterm | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).SelfControlFailureLongterm))
                names{3}     = 'Self Control Failure | Longterm';
                onsets{3}    = CDM_FU3.CDM_Onsets(j).SelfControlFailureLongterm;
                durations{3} = zeros(numel(CDM_FU3.CDM_Onsets(j).SelfControlFailureLongterm),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: SelfControlFailureLongterm | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).SelfControlFailureShortterm))
                names{4}     = 'Self Control Failure | Shortterm';
                onsets{4}    = CDM_FU3.CDM_Onsets(j).SelfControlFailureShortterm;
                durations{4} = zeros(numel(CDM_FU3.CDM_Onsets(j).SelfControlFailureShortterm),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: SelfControlFailureShortterm | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).NoConflictAllPositiveYes))
                names{5}     = 'No Conflict | All Positive | Accepted';
                onsets{5}    = CDM_FU3.CDM_Onsets(j).NoConflictAllPositiveYes;
                durations{5} = zeros(numel(CDM_FU3.CDM_Onsets(j).NoConflictAllPositiveYes),1);
                orth{5}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: NoConflictAllPositiveYes | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).NoConflictAllPositiveNo))
                names{6}     = 'No Conflict | All Positive | Declined';
                onsets{6}    = CDM_FU3.CDM_Onsets(j).NoConflictAllPositiveNo;
                durations{6} = zeros(numel(CDM_FU3.CDM_Onsets(j).NoConflictAllPositiveNo),1);
                orth{6}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: NoConflictAllPositiveNo | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).NoConflictAllNegativeYes))
                names{7}     = 'No Conflict | All Negative | Accepted';
                onsets{7}    = CDM_FU3.CDM_Onsets(j).NoConflictAllNegativeYes;
                durations{7} = zeros(numel(CDM_FU3.CDM_Onsets(j).NoConflictAllNegativeYes),1);
                orth{7}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: NoConflictAllNegativeYes | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU3.CDM_Onsets(j).NoConflictAllNegativeNo))
                names{8}     = 'No Conflict | All Negative | Declined';
                onsets{8}    = CDM_FU3.CDM_Onsets(j).NoConflictAllNegativeNo;
                durations{8} = zeros(numel(CDM_FU3.CDM_Onsets(j).NoConflictAllNegativeNo),1);
                orth{8}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU3 | No onset of: NoConflictAllNegativeNo | Subject: %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                n = n+1;
            end              

            if exist('onsets','var')
                %% Remove empty cells (if they exist)
                names     = names(~cellfun('isempty',names));
                onsets    = onsets(~cellfun('isempty',onsets));
                durations = durations(~cellfun('isempty',durations));
                orth      = orth(~cellfun('isempty',orth));
                %% save onset files
                save([outpath num2str(CDM_FU3.CDM_Onsets(j).SubjectID) '_CDM_FU3_mk2021' '.mat'], 'names','onsets', 'durations','orth');
                %% purge variables
                clear names onsets durations orth
            else
                MissingDataSet{m,1} = sprintf('Due to missing data, cannot create onset file | _CDM_FU3_mk2021 | %s', num2str(CDM_FU3.CDM_Onsets(j).SubjectID));
                m = m+1;
            end

        end

    elseif strcmp(sessTasks{i},'CDM_FU6')

        IDs = [CDM_FU6.CDM_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Complex Decision Making task (Sixth follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %%
            if isempty(CDM_FU6.CDM_Onsets(j).SubjectID)
                continue % skip iteration if data set for that is missing
            end

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2020 | 1st GLM

            %%
            names{1}         = 'Stimulus Decision';
            onsets{1}        = CDM_FU6.CDM_Onsets(j).StimulusDecision;
            durations{1}     = CDM_FU6.CDM_Onsets(j).RT_StimulusDecision;
            pmod(1).name{1}  = 'DecisionValue';
            pmod(1).param{1} = CDM_FU6.CDM_Onsets(j).StimulusDecisionParametricReversed;
            pmod(1).poly{1}  = 1;
            orth{1}          = true;

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).MissedTrials))
                names{2}     = 'Missed Trials';
                onsets{2}    = CDM_FU6.CDM_Onsets(j).MissedTrials;
                durations{2} = zeros(numel(CDM_FU6.CDM_Onsets(j).MissedTrials),1);
                orth{2}      = true;
            end

            %% save onset files
            save([outpath num2str(CDM_FU6.CDM_Onsets(j).SubjectID) '_CDM_FU6_mk2020_GLM01' '.mat'], 'names','onsets', 'durations','pmod','orth');
            %% purge variables
            clear pmod

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2020 | 2nd GLM
            pmod(1).name{1}  = 'StimulusRatingReversed.ShortTerm';
            pmod(1).param{1} = CDM_FU6.CDM_Onsets(j).StimulusDecisionParametricRatingsReversed(:,1);
            pmod(1).poly{1}  = 1;
            pmod(1).name{2}  = 'StimulusRatingReversed.LongTerm';
            pmod(1).param{2} = CDM_FU6.CDM_Onsets(j).StimulusDecisionParametricRatingsReversed(:,2);
            pmod(1).poly{2}  = 1;

            %% save onset files
            save([outpath num2str(CDM_FU6.CDM_Onsets(j).SubjectID) '_CDM_FU6_mk2020_GLM02' '.mat'], 'names','onsets', 'durations','pmod','orth');
            %% purge variables
            clear names onsets durations pmod orth

            %% Create onset files for Complex Decision Making task | akin to Kroenke et al., 2021

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).SelfControlSuccessLongterm))
                names{1}     = 'Self Control Success | Longterm';
                onsets{1}    = CDM_FU6.CDM_Onsets(j).SelfControlSuccessLongterm;
                durations{1} = zeros(numel(CDM_FU6.CDM_Onsets(j).SelfControlSuccessLongterm),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: SelfControlSuccessLongterm | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).SelfControlSuccessShortterm))
                names{2}     = 'Self Control Success | Shortterm';
                onsets{2}    = CDM_FU6.CDM_Onsets(j).SelfControlSuccessShortterm;
                durations{2} = zeros(numel(CDM_FU6.CDM_Onsets(j).SelfControlSuccessShortterm),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: SelfControlSuccessShortterm | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).SelfControlFailureLongterm))
                names{3}     = 'Self Control Failure | Longterm';
                onsets{3}    = CDM_FU6.CDM_Onsets(j).SelfControlFailureLongterm;
                durations{3} = zeros(numel(CDM_FU6.CDM_Onsets(j).SelfControlFailureLongterm),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: SelfControlFailureLongterm | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).SelfControlFailureShortterm))
                names{4}     = 'Self Control Failure | Shortterm';
                onsets{4}    = CDM_FU6.CDM_Onsets(j).SelfControlFailureShortterm;
                durations{4} = zeros(numel(CDM_FU6.CDM_Onsets(j).SelfControlFailureShortterm),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: SelfControlFailureShortterm | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).NoConflictAllPositiveYes))
                names{5}     = 'No Conflict | All Positive | Accepted';
                onsets{5}    = CDM_FU6.CDM_Onsets(j).NoConflictAllPositiveYes;
                durations{5} = zeros(numel(CDM_FU6.CDM_Onsets(j).NoConflictAllPositiveYes),1);
                orth{5}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: NoConflictAllPositiveYes | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).NoConflictAllPositiveNo))
                names{6}     = 'No Conflict | All Positive | Declined';
                onsets{6}    = CDM_FU6.CDM_Onsets(j).NoConflictAllPositiveNo;
                durations{6} = zeros(numel(CDM_FU6.CDM_Onsets(j).NoConflictAllPositiveNo),1);
                orth{6}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: NoConflictAllPositiveNo | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).NoConflictAllNegativeYes))
                names{7}     = 'No Conflict | All Negative | Accepted';
                onsets{7}    = CDM_FU6.CDM_Onsets(j).NoConflictAllNegativeYes;
                durations{7} = zeros(numel(CDM_FU6.CDM_Onsets(j).NoConflictAllNegativeYes),1);
                orth{7}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: NoConflictAllNegativeYes | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(CDM_FU6.CDM_Onsets(j).NoConflictAllNegativeNo))
                names{8}     = 'No Conflict | All Negative | Declined';
                onsets{8}    = CDM_FU6.CDM_Onsets(j).NoConflictAllNegativeNo;
                durations{8} = zeros(numel(CDM_FU6.CDM_Onsets(j).NoConflictAllNegativeNo),1);
                orth{8}      = true;
            else
                MissingMessage{n,1} = sprintf('CDM | FU6 | No onset of: NoConflictAllNegativeNo | Subject: %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                n = n+1;
            end

           
            if exist('onsets','var')
                %% Remove empty cells (if they exist)
                names     = names(~cellfun('isempty',names));
                onsets    = onsets(~cellfun('isempty',onsets));
                durations = durations(~cellfun('isempty',durations));
                orth      = orth(~cellfun('isempty',orth));
                %% save onset files
                save([outpath num2str(CDM_FU6.CDM_Onsets(j).SubjectID) '_CDM_FU6_mk2021' '.mat'], 'names','onsets', 'durations','orth');
                %% purge variables
                clear names onsets durations orth
            else
                MissingDataSet{m,1} = sprintf('Due to missing data, cannot create onset file | _CDM_FU6_mk2021 | %s', num2str(CDM_FU6.CDM_Onsets(j).SubjectID));
                m = m+1;
            end

        end

    elseif strcmp(sessTasks{i},'GoNogo_BL')

        IDs = [GoNogo_BL.GoNogo_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Go/Nogo task (Baseline session)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %% Create onset files for the Go/Nogo task |

            %%
            if ~any(isnan(GoNogo_BL.GoNogo_Onsets(j).GoCorrect))
                names{1}     = 'Go Correct';
                onsets{1}    = GoNogo_BL.GoNogo_Onsets(j).GoCorrect;
                durations{1} = zeros(numel(GoNogo_BL.GoNogo_Onsets(j).GoCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | BL | No onset of: GoCorrect | Subject: %s', num2str(GoNogo_BL.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_BL.GoNogo_Onsets(j).NogoCorrect))
                names{2}     = 'Nogo Correct';
                onsets{2}    = GoNogo_BL.GoNogo_Onsets(j).NogoCorrect;
                durations{2} = zeros(numel(GoNogo_BL.GoNogo_Onsets(j).NogoCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | BL | No onset of: NogoCorrect | Subject: %s', num2str(GoNogo_BL.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_BL.GoNogo_Onsets(j).NogoError))
                names{3}     = 'Nogo Error';
                onsets{3}    = GoNogo_BL.GoNogo_Onsets(j).NogoError;
                durations{3} = zeros(numel(GoNogo_BL.GoNogo_Onsets(j).NogoError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | BL | No onset of: NogoError | Subject: %s', num2str(GoNogo_BL.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_BL.GoNogo_Onsets(j).GoError))
                names{4}     = 'Go Error';
                onsets{4}    = GoNogo_BL.GoNogo_Onsets(j).GoError;
                durations{4} = zeros(numel(GoNogo_BL.GoNogo_Onsets(j).GoError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | BL | No onset of: GoError | Subject: %s', num2str(GoNogo_BL.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(GoNogo_BL.GoNogo_Onsets(j).SubjectID) '_GNG_BL' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'GoNogo_FU3')

        IDs = [GoNogo_FU3.GoNogo_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Go/Nogo task (Third follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %% Create onset files for the Go/Nogo task |

            %%
            if ~any(isnan(GoNogo_FU3.GoNogo_Onsets(j).GoCorrect))
                names{1}     = 'Go Correct';
                onsets{1}    = GoNogo_FU3.GoNogo_Onsets(j).GoCorrect;
                durations{1} = zeros(numel(GoNogo_FU3.GoNogo_Onsets(j).GoCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU3 | No onset of: GoCorrect | Subject: %s', num2str(GoNogo_FU3.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_FU3.GoNogo_Onsets(j).NogoCorrect))
                names{2}     = 'Nogo Correct';
                onsets{2}    = GoNogo_FU3.GoNogo_Onsets(j).NogoCorrect;
                durations{2} = zeros(numel(GoNogo_FU3.GoNogo_Onsets(j).NogoCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU3 | No onset of: NogoCorrect | Subject: %s', num2str(GoNogo_FU3.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_FU3.GoNogo_Onsets(j).NogoError))
                names{3}     = 'Nogo Error';
                onsets{3}    = GoNogo_FU3.GoNogo_Onsets(j).NogoError;
                durations{3} = zeros(numel(GoNogo_FU3.GoNogo_Onsets(j).NogoError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU3 | No onset of: NogoError | Subject: %s', num2str(GoNogo_FU3.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_FU3.GoNogo_Onsets(j).GoError))
                names{4}     = 'Go Error';
                onsets{4}    = GoNogo_FU3.GoNogo_Onsets(j).GoError;
                durations{4} = zeros(numel(GoNogo_FU3.GoNogo_Onsets(j).GoError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU3 | No onset of: GoError | Subject: %s', num2str(GoNogo_FU3.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(GoNogo_FU3.GoNogo_Onsets(j).SubjectID) '_GNG_FU3' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'GoNogo_FU6')

        IDs = [GoNogo_FU6.GoNogo_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Go/Nogo task (Sixth follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %% Create onset files for the Go/Nogo task

            %%
            if ~any(isnan(GoNogo_FU6.GoNogo_Onsets(j).GoCorrect))
                names{1}     = 'Go Correct';
                onsets{1}    = GoNogo_FU6.GoNogo_Onsets(j).GoCorrect;
                durations{1} = zeros(numel(GoNogo_FU6.GoNogo_Onsets(j).GoCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU6 | No onset of: GoCorrect | Subject: %s', num2str(GoNogo_FU6.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_FU6.GoNogo_Onsets(j).NogoCorrect))
                names{2}     = 'Nogo Correct';
                onsets{2}    = GoNogo_FU6.GoNogo_Onsets(j).NogoCorrect;
                durations{2} = zeros(numel(GoNogo_FU6.GoNogo_Onsets(j).NogoCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU6 | No onset of: NogoCorrect | Subject: %s', num2str(GoNogo_FU6.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_FU6.GoNogo_Onsets(j).NogoError))
                names{3}     = 'Nogo Error';
                onsets{3}    = GoNogo_FU6.GoNogo_Onsets(j).NogoError;
                durations{3} = zeros(numel(GoNogo_FU6.GoNogo_Onsets(j).NogoError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU6 | No onset of: NogoError | Subject: %s', num2str(GoNogo_FU6.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(GoNogo_FU6.GoNogo_Onsets(j).GoError))
                names{4}     = 'Go Error';
                onsets{4}    = GoNogo_FU6.GoNogo_Onsets(j).GoError;
                durations{4} = zeros(numel(GoNogo_FU6.GoNogo_Onsets(j).GoError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('GNG | FU6 | No onset of: GoError | Subject: %s', num2str(GoNogo_FU6.GoNogo_Onsets(j).SubjectID));
                n = n+1;
            end


            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(GoNogo_FU6.GoNogo_Onsets(j).SubjectID) '_GNG_FU6' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'Stroop_BL')

        IDs = [Stroop_BL.Stroop_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Counting Stroop task (Baseline session)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %% Create onset files for the Counting Stroop task | akin to Kroenke et al., 2015

            %%
            if ~any(isnan(Stroop_BL.Stroop_Onsets(j).ConCorrect))
                names{1}     = 'Congruent Trials (correct)';
                onsets{1}    = Stroop_BL.Stroop_Onsets(j).ConCorrect;
                durations{1} = zeros(numel(Stroop_BL.Stroop_Onsets(j).ConCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | BL | No onset of: ConCorrect | Subject: %s', num2str(Stroop_BL.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_BL.Stroop_Onsets(j).IncCorrect))
                names{2}     = 'Incongruent Trials (correct)';
                onsets{2}    = Stroop_BL.Stroop_Onsets(j).IncCorrect;
                durations{2} = zeros(numel(Stroop_BL.Stroop_Onsets(j).IncCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | BL | No onset of: IncCorrect | Subject: %s', num2str(Stroop_BL.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_BL.Stroop_Onsets(j).IncError))
                names{3}     = 'Incongruent errors';
                onsets{3}    = Stroop_BL.Stroop_Onsets(j).IncError;
                durations{3} = zeros(numel(Stroop_BL.Stroop_Onsets(j).IncError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | BL | No onset of: IncError | Subject: %s', num2str(Stroop_BL.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_BL.Stroop_Onsets(j).ConError))
                names{4}     = 'Congruent errors';
                onsets{4}    = Stroop_BL.Stroop_Onsets(j).ConError;
                durations{4} = zeros(numel(Stroop_BL.Stroop_Onsets(j).ConError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | BL | No onset of: ConError | Subject: %s', num2str(Stroop_BL.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(Stroop_BL.Stroop_Onsets(j).SubjectID) '_STR_BL_mk2015' '.mat'], 'names','onsets', 'durations','orth');


            %% Create onset files for the Counting Stroop task | akin to Kroenke et al., 2018

            %%
            if ~any(isnan(Stroop_BL.Stroop_Onsets(j).PostIncError))
                names{5}     = 'Post-incongruent errors';
                onsets{5}    = Stroop_BL.Stroop_Onsets(j).PostIncError;
                durations{5} = zeros(numel(Stroop_BL.Stroop_Onsets(j).PostIncError),1);
                orth{5}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | BL | No onset of: PostIncError | Subject: %s', num2str(Stroop_BL.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_BL.Stroop_Onsets(j).PostIncCorrect))
                names{6}     = 'Post-incongruent correct trials';
                onsets{6}    = Stroop_BL.Stroop_Onsets(j).PostIncCorrect;
                durations{6} = zeros(numel(Stroop_BL.Stroop_Onsets(j).PostIncCorrect),1);
                orth{6}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | BL | No onset of: PostIncCorrect | Subject: %s', num2str(Stroop_BL.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end


            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(Stroop_BL.Stroop_Onsets(j).SubjectID) '_STR_BL_mk2018' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'Stroop_FU3')

        IDs = [Stroop_FU3.Stroop_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Counting Stroop task (Third follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %% Create onset files for the Counting Stroop task | akin to Kroenke et al., 2015

            %%
            if ~any(isnan(Stroop_FU3.Stroop_Onsets(j).ConCorrect))
                names{1}     = 'Congruent Trials (correct)';
                onsets{1}    = Stroop_FU3.Stroop_Onsets(j).ConCorrect;
                durations{1} = zeros(numel(Stroop_FU3.Stroop_Onsets(j).ConCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU3 | No onset of: ConCorrect | Subject: %s', num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU3.Stroop_Onsets(j).IncCorrect))
                names{2}     = 'Incongruent Trials (correct)';
                onsets{2}    = Stroop_FU3.Stroop_Onsets(j).IncCorrect;
                durations{2} = zeros(numel(Stroop_FU3.Stroop_Onsets(j).IncCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU3 | No onset of: IncCorrect | Subject: %s', num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU3.Stroop_Onsets(j).IncError))
                names{3}     = 'Incongruent errors';
                onsets{3}    = Stroop_FU3.Stroop_Onsets(j).IncError;
                durations{3} = zeros(numel(Stroop_FU3.Stroop_Onsets(j).IncError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | BL | No onset of: IncError | Subject: %s', num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU3.Stroop_Onsets(j).ConError))
                names{4}     = 'Congruent errors';
                onsets{4}    = Stroop_FU3.Stroop_Onsets(j).ConError;
                durations{4} = zeros(numel(Stroop_FU3.Stroop_Onsets(j).ConError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU3 | No onset of: ConError | Subject: %s', num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID) '_STR_FU3_mk2015' '.mat'], 'names','onsets', 'durations','orth');


            %% Create onset files for the Counting Stroop task | akin to Kroenke et al., 2018

            %%
            if ~any(isnan(Stroop_FU3.Stroop_Onsets(j).PostIncError))
                names{5}     = 'Post-incongruent errors';
                onsets{5}    = Stroop_FU3.Stroop_Onsets(j).PostIncError;
                durations{5} = zeros(numel(Stroop_FU3.Stroop_Onsets(j).PostIncError),1);
                orth{5}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU3 | No onset of: PostIncError | Subject: %s', num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU3.Stroop_Onsets(j).PostIncCorrect))
                names{6}     = 'Post-incongruent correct trials';
                onsets{6}    = Stroop_FU3.Stroop_Onsets(j).PostIncCorrect;
                durations{6} = zeros(numel(Stroop_FU3.Stroop_Onsets(j).PostIncCorrect),1);
                orth{6}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU3 | No onset of: PostIncCorrect | Subject: %s', num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end


            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(Stroop_FU3.Stroop_Onsets(j).SubjectID) '_STR_FU3_mk2018' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'Stroop_FU6')

        IDs = [Stroop_FU6.Stroop_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Counting Stroop task (Sixth follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);

            %% Create onset files for the Counting Stroop task | akin to Kroenke et al., 2015

            %%
            if ~any(isnan(Stroop_FU6.Stroop_Onsets(j).ConCorrect))
                names{1}     = 'Congruent Trials (correct)';
                onsets{1}    = Stroop_FU6.Stroop_Onsets(j).ConCorrect;
                durations{1} = zeros(numel(Stroop_FU6.Stroop_Onsets(j).ConCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU6 | No onset of: ConCorrect | Subject: %s', num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU6.Stroop_Onsets(j).IncCorrect))
                names{2}     = 'Incongruent Trials (correct)';
                onsets{2}    = Stroop_FU6.Stroop_Onsets(j).IncCorrect;
                durations{2} = zeros(numel(Stroop_FU6.Stroop_Onsets(j).IncCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU6 | No onset of: IncCorrect | Subject: %s', num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU6.Stroop_Onsets(j).IncError))
                names{3}     = 'Incongruent errors';
                onsets{3}    = Stroop_FU6.Stroop_Onsets(j).IncError;
                durations{3} = zeros(numel(Stroop_FU6.Stroop_Onsets(j).IncError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU6 | No onset of: IncError | Subject: %s', num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU6.Stroop_Onsets(j).ConError))
                names{4}     = 'Congruent errors';
                onsets{4}    = Stroop_FU6.Stroop_Onsets(j).ConError;
                durations{4} = zeros(numel(Stroop_FU6.Stroop_Onsets(j).ConError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU6 | No onset of: ConError | Subject: %s', num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID) '_STR_FU6_mk2015' '.mat'], 'names','onsets', 'durations','orth');


            %% Create onset files for the Counting Stroop task | akin to Kroenke et al., 2018

            %%
            if ~any(isnan(Stroop_FU6.Stroop_Onsets(j).PostIncError))
                names{5}     = 'Post-incongruent errors';
                onsets{5}    = Stroop_FU6.Stroop_Onsets(j).PostIncError;
                durations{5} = zeros(numel(Stroop_FU6.Stroop_Onsets(j).PostIncError),1);
                orth{5}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU6 | No onset of: PostIncError | Subject: %s', num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(Stroop_FU6.Stroop_Onsets(j).PostIncCorrect))
                names{6}     = 'Post-incongruent correct trials';
                onsets{6}    = Stroop_FU6.Stroop_Onsets(j).PostIncCorrect;
                durations{6} = zeros(numel(Stroop_FU6.Stroop_Onsets(j).PostIncCorrect),1);
                orth{6}      = true;
            else
                MissingMessage{n,1} = sprintf('STR | FU6 | No onset of: PostIncCorrect | Subject: %s', num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID));
                n = n+1;
            end


            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(Stroop_FU6.Stroop_Onsets(j).SubjectID) '_STR_FU6_mk2018' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'StopSignal_BL')

        IDs = [StopSignal_BL.StopSignal_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Stop-Signal task (Baseline session)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);


            %% Create onset files for the Stop-Signal task

            %%
            if ~any(isnan(StopSignal_BL.StopSignal_Onsets(j).GoCorrect))
                names{1}     = 'Go Trials (correct)';
                onsets{1}    = StopSignal_BL.StopSignal_Onsets(j).GoCorrect;
                durations{1} = zeros(numel(StopSignal_BL.StopSignal_Onsets(j).GoCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | BL | No onset of: GoCorrect | Subject: %s', num2str(StopSignal_BL.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_BL.StopSignal_Onsets(j).StopCorrect))
                names{2}     = 'Stop Trials (correct)';
                onsets{2}    = StopSignal_BL.StopSignal_Onsets(j).StopCorrect;
                durations{2} = zeros(numel(StopSignal_BL.StopSignal_Onsets(j).StopCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | BL | No onset of: StopCorrect | Subject: %s', num2str(StopSignal_BL.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_BL.StopSignal_Onsets(j).StopError))
                names{3}     = 'Stop errors';
                onsets{3}    = StopSignal_BL.StopSignal_Onsets(j).StopError;
                durations{3} = zeros(numel(StopSignal_BL.StopSignal_Onsets(j).StopError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | BL | No onset of: StopError | Subject: %s', num2str(StopSignal_BL.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_BL.StopSignal_Onsets(j).GoError))
                names{4}     = 'Go errors';
                onsets{4}    = StopSignal_BL.StopSignal_Onsets(j).GoError;
                durations{4} = zeros(numel(StopSignal_BL.StopSignal_Onsets(j).GoError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | BL | No onset of: GoError | Subject: %s', num2str(StopSignal_BL.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(StopSignal_BL.StopSignal_Onsets(j).SubjectID) '_SST_BL' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'StopSignal_FU3')

        IDs = [StopSignal_FU3.StopSignal_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Stop-Signal task (Third follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);


            %% Create onset files for the Stop-Signal task

            %%
            if ~any(isnan(StopSignal_FU3.StopSignal_Onsets(j).GoCorrect))
                names{1}     = 'Go Trials (correct)';
                onsets{1}    = StopSignal_FU3.StopSignal_Onsets(j).GoCorrect;
                durations{1} = zeros(numel(StopSignal_FU3.StopSignal_Onsets(j).GoCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU3 | No onset of: GoCorrect | Subject: %s', num2str(StopSignal_FU3.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_FU3.StopSignal_Onsets(j).StopCorrect))
                names{2}     = 'Stop Trials (correct)';
                onsets{2}    = StopSignal_FU3.StopSignal_Onsets(j).StopCorrect;
                durations{2} = zeros(numel(StopSignal_FU3.StopSignal_Onsets(j).StopCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU3 | No onset of: StopCorrect | Subject: %s', num2str(StopSignal_FU3.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_FU3.StopSignal_Onsets(j).StopError))
                names{3}     = 'Stop errors';
                onsets{3}    = StopSignal_FU3.StopSignal_Onsets(j).StopError;
                durations{3} = zeros(numel(StopSignal_FU3.StopSignal_Onsets(j).StopError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU3 | No onset of: StopError | Subject: %s', num2str(StopSignal_FU3.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_FU3.StopSignal_Onsets(j).GoError))
                names{4}     = 'Go errors';
                onsets{4}    = StopSignal_FU3.StopSignal_Onsets(j).GoError;
                durations{4} = zeros(numel(StopSignal_FU3.StopSignal_Onsets(j).GoError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU3 | No onset of: GoError | Subject: %s', num2str(StopSignal_FU3.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end


            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(StopSignal_FU3.StopSignal_Onsets(j).SubjectID) '_SST_FU3' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end

    elseif strcmp(sessTasks{i},'StopSignal_FU6')

        IDs = [StopSignal_FU6.StopSignal_Onsets.SubjectID]';
        %%
        fprintf(1,'\nProcessing stimulus onset files of Stop-Signal task (Sixth follow-up)\n')
        fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
        fprintf(1,'%d datasets in total\n',numel(IDs))
        fprintf(1,'Processing subject #: ');

        for j = 1:numel(IDs)
            %%
            if j>1
                for b=0:log10(j-1)
                    fprintf(1,'\b');  % Deleting last displayed number
                end
            end
            fprintf(1,'%d',j);


            %% Create onset files for the Stop-Signal task

            %%
            if ~any(isnan(StopSignal_FU6.StopSignal_Onsets(j).GoCorrect))
                names{1}     = 'Go Trials (correct)';
                onsets{1}    = StopSignal_FU6.StopSignal_Onsets(j).GoCorrect;
                durations{1} = zeros(numel(StopSignal_FU6.StopSignal_Onsets(j).GoCorrect),1);
                orth{1}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU6 | No onset of: GoCorrect | Subject: %s', num2str(StopSignal_FU6.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_FU6.StopSignal_Onsets(j).StopCorrect))
                names{2}     = 'Stop Trials (correct)';
                onsets{2}    = StopSignal_FU6.StopSignal_Onsets(j).StopCorrect;
                durations{2} = zeros(numel(StopSignal_FU6.StopSignal_Onsets(j).StopCorrect),1);
                orth{2}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU6 | No onset of: StopCorrect | Subject: %s', num2str(StopSignal_FU6.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_FU6.StopSignal_Onsets(j).StopError))
                names{3}     = 'Stop errors';
                onsets{3}    = StopSignal_FU6.StopSignal_Onsets(j).StopError;
                durations{3} = zeros(numel(StopSignal_FU6.StopSignal_Onsets(j).StopError),1);
                orth{3}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU6 | No onset of: StopError | Subject: %s', num2str(StopSignal_FU6.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %%
            if ~any(isnan(StopSignal_FU6.StopSignal_Onsets(j).GoError))
                names{4}     = 'Go errors';
                onsets{4}    = StopSignal_FU6.StopSignal_Onsets(j).GoError;
                durations{4} = zeros(numel(StopSignal_FU6.StopSignal_Onsets(j).GoError),1);
                orth{4}      = true;
            else
                MissingMessage{n,1} = sprintf('SST | FU6 | No onset of: GoError | Subject: %s', num2str(StopSignal_FU6.StopSignal_Onsets(j).SubjectID));
                n = n+1;
            end

            %% Remove empty cells (if they exist)
            names     = names(~cellfun('isempty',names));
            onsets    = onsets(~cellfun('isempty',onsets));
            durations = durations(~cellfun('isempty',durations));
            orth      = orth(~cellfun('isempty',orth));
            %% save onset files
            save([outpath num2str(StopSignal_FU6.StopSignal_Onsets(j).SubjectID) '_SST_FU6' '.mat'], 'names','onsets', 'durations','orth');
            %% purge variables
            clear names onsets durations orth

        end
    end
end

%%
for k = 1:numel(MissingMessage)
    tmp                  = strsplit(MissingMessage{k,1},'|');
    Task{k,1}            = deblank(tmp(1));
    Session{k,1}         = deblank(tmp(2));
    MissingVariable{k,1} = deblank(tmp{3}(15:end));
    Subject{k,1}         = tmp{4}(11:end);
    clear tmp
end

%%
MissingTable = table(Task,Session,MissingVariable,Subject);
MissingTable.Properties.VariableNames = {'Task','Session','MissingVariable','Subject'};
writetable(MissingTable,[outpath 'MissingVariablesTable.xlsx']);
fprintf('\n')
disp(MissingTable);
fprintf('\n')
%%
