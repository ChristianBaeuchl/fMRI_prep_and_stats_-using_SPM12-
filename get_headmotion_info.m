function get_headmotion_info
%%get_headmotion_info
% Creates an overview of translational and rotational head motion outliers
% based on predefined criteria and plots the progression of head motion 
%
% Christian Baeuchl | 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear variables, close windows and clear the command line
clc, clear, close
%%


% define specifics of data to be examined
task = {'task name','task abbreviation'};
sess = 'session indicator';

% define paths and change path
path       = ['path to your preprocessed data' sess 'uses session information'];
returnpath = 'path location to where the folder with head motion plots should be created';
%cd(path)

% load data list
list = sort(importdata(['.\Participant_lists\List_' task{2} '_' sess '.txt']));

% outlier criteria
driftcrit = 3;
spikecrit = 0.5;

% define in- and output folders
functfold = [task{1} '_EPI'];
type      = exist([returnpath filesep 'HeadMotionPlots' filesep task{1} '_' sess(4:end) ],'dir');
if type ~= 7
    mkdir([returnpath filesep 'HeadMotionPlots' filesep task{1} '_' sess(4:end)]);
end

%% call function to calculate Volterra expansion of the realignment parameters
volt_exp_multsub(task,path,list)
%%

% which realignment parameter file should be used?
rpfile = 'Vxrp*.txt'; %rp*.txt

%%
DrifttransoutX = nan(length(list),1);
DrifttransoutY = nan(length(list),1);
DrifttransoutZ = nan(length(list),1);
DriftrotoutX   = cellstr(char((nan(length(list),1))));
DriftrotoutY   = cellstr(char((nan(length(list),1))));
DriftrotoutZ   = cellstr(char((nan(length(list),1))));

SpiketransoutX = nan(length(list),1);
SpiketransoutY = nan(length(list),1);
SpiketransoutZ = nan(length(list),1);
SpikerotoutX   = nan(length(list),1);
SpikerotoutY   = nan(length(list),1);
SpikerotoutZ   = nan(length(list),1);

SpiketransoutXPerc = nan(length(list),1);
SpiketransoutYPerc = nan(length(list),1);
SpiketransoutZPerc = nan(length(list),1);
SpikerotoutXPerc   = nan(length(list),1);
SpikerotoutYPerc   = nan(length(list),1);
SpikerotoutZPerc   = nan(length(list),1);

NumSamples         = nan(length(list),1);

%%
fprintf(1,'\nChecking and plotting (head) motion parameters\n')
fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  \n');
fprintf(1,'\n');
fprintf(1,'%d data sets in total\n',size(list,1))
fprintf(1,'Processing data set #: ' );

%% loop
for i = 1:length(list)

    %%
    if i>1
        for b=0:log10(i-1)
            fprintf(1,'\b');  % Deleting last displayed number
        end
    end
    fprintf(1,'%d',i);

    cd([path num2str(list(i)) filesep functfold])
    rpname = ls(rpfile);
    rp     = load(rpname);
    NumSamples(i) = size(rp,1);

    trans = find(max(abs(rp(:,1:3))) > driftcrit);
    if ~isempty(trans(trans==1)) % x - translation drift outlier
        DrifttransoutX(i) = 1;
    end
    if ~isempty(trans(trans==2)) % y - translation drift outlier
        DrifttransoutY(i) = 1;
    end
    if ~isempty(trans(trans==3)) % z - translation drift outlier
        DrifttransoutZ(i) = 1;
    end

    rot   = find(max(abs(rp(:,4:6)*(180/pi))) > driftcrit);
    if ~isempty(rot(rot==1)) % x - rotation drift outlier
        DriftrotoutX(i) = {'pitch'};
    end
    if ~isempty(rot(rot==2)) % y - rotation drift outlier
        DriftrotoutY(i) = {'roll'};
    end
    if ~isempty(rot(rot==3)) % z - rotation drift outlier
        DriftrotoutZ(i) = {'yaw'};
    end

    [~,spiket] = find(abs(rp(:,7:9)) > spikecrit);
    numspiketx = numel(spiket(spiket==1));
    numspikety = numel(spiket(spiket==2));
    numspiketz = numel(spiket(spiket==3));

    if ~isempty(numspiketx) % x - translation drift outlier
        SpiketransoutX(i) = numspiketx;
    end
    if ~isempty(numspikety) % y - translation drift outlier
        SpiketransoutY(i) = numspikety;
    end
    if ~isempty(numspiketz) % z - translation drift outlier
        SpiketransoutZ(i) = numspiketz;
    end

    [~,spiker] = find((abs(rp(:,10:12))*(180/pi)) > spikecrit);
    numspikerx = numel(spiker(spiker==1));
    numspikery = numel(spiker(spiker==2));
    numspikerz = numel(spiker(spiker==3));

    if ~isempty(numspikerx) % x - rotation drift outlier
        SpikerotoutX(i) = numspikerx;
    end
    if ~isempty(numspikety) % y - rotation drift outlier
        SpikerotoutY(i) = numspikery;
    end
    if ~isempty(numspiketz) % z - rotation drift outlier
        SpikerotoutZ(i) = numspikerz;
    end

    %%
    SpiketransoutXPerc(i) = (SpiketransoutX(i)*100)/NumSamples(i);
    SpiketransoutYPerc(i) = (SpiketransoutY(i)*100)/NumSamples(i);
    SpiketransoutZPerc(i) = (SpiketransoutZ(i)*100)/NumSamples(i);
    SpikerotoutXPerc(i)   = (SpikerotoutX(i)*100)/NumSamples(i);
    SpikerotoutYPerc(i)   = (SpikerotoutY(i)*100)/NumSamples(i);
    SpikerotoutZPerc(i)   = (SpikerotoutZ(i)*100)/NumSamples(i);


    %% plotting

    figure;
    set(gcf,'color','w')
    subplot(2,2,1)
    plot(rp(:,1),'-b')
    hold on
    grid on
    plot(rp(:,2),'-g')
    plot(rp(:,3),'-r')
    title([num2str(list(i)) '- Translation']);

    xlabel('time (scans)')
    ylabel('amplitude (milimeter)')
    legend('x','y','z')
    axis([1 length(rp) min(min(rp(:,1:3)))-1 max(max(rp(:,1:3)))+1])
    subplot(2,2,3)
    plot(rp(:,7),'-b')
    hold on
    grid on
    plot(rp(:,8),'-g')
    plot(rp(:,9),'-r')
    plot(spikecrit*ones(1,length(rp)),'m');  % line for outliers
    plot(-spikecrit*ones(1,length(rp)),'m'); % line for outliers
    title([num2str(list(i)) '- Translation (relative displacement)']);
    xlabel('time (scans)')
    ylabel('amplitude (milimeters)')
    legend('x','y','z')
    axis([1 length(rp) min(min(rp(:,7:9)))-1 max(max(rp(:,7:9)))+1])
    hold off

    subplot(2,2,2)
    plot(rp(:,4)*(180/pi),'-b')
    hold on
    grid on
    plot(rp(:,5)*(180/pi),'-g')
    plot(rp(:,6)*(180/pi),'-r')
    title([num2str(list(i)) '- Rotation']);
    xlabel('time (scans)')
    ylabel('degree')
    legend('pitch','roll','yaw')
    axis([1 length(rp) min(min(rp(:,4:6)*(180/pi)))-1 max(max(rp(:,4:6)*(180/pi)))+1])
    subplot(2,2,4)
    plot(rp(:,10)*(180/pi),'-b')
    hold on
    grid on
    plot(rp(:,11)*(180/pi),'-g')
    plot(rp(:,12)*(180/pi),'-r')
    plot(spikecrit*ones(1,length(rp)),'m'); % line for outliers
    plot(-spikecrit*ones(1,length(rp)),'m'); % line for outliers
    title([num2str(list(i)) '- Rotation (relative displacement)']);
    xlabel('time (scans)')
    ylabel('degree')
    legend('pitch','roll','yaw')
    axis([1 length(rp) min(min(rp(:,10:12)*(180/pi)))-1 max(max(rp(:,10:12)*(180/pi)))+1])
    set(gcf,'PaperOrientation','landscape');
    hold off

    cd([returnpath 'HeadMotionPlots' filesep]);
    print(num2str(list(i)),'-fillpage','-dpdf');
    close(gcf)

    clear rp
end
cd(returnpath)
fprintf('\n')


%%
MotTable = table(list,DrifttransoutX,DrifttransoutY,DrifttransoutZ,DriftrotoutX,DriftrotoutY,DriftrotoutZ,...
    SpiketransoutX,SpiketransoutY,SpiketransoutZ,SpikerotoutX,SpikerotoutY,SpikerotoutZ,NumSamples,...
    SpiketransoutXPerc,SpiketransoutYPerc,SpiketransoutZPerc,SpikerotoutXPerc,SpikerotoutYPerc,SpikerotoutZPerc);
MotTable.Properties.VariableNames = {'ID','TranslationDriftX','TranslationDriftY','TranslationDriftZ','RotationDriftX','RotationDriftY',...
    'RotationDriftZ','TranslationSpikesX','TranslationSpikesY','TranslationSpikesZ','RotationSpikesX','RotationSpikesY','RotationSpikesZ',...
    'NumOfSamples','TranslationSpikesXPerc','TranslationSpikesYPerc','TranslationSpikesZPerc','RotationSpikesXPerc','RotationSpikesYPerc',...
    'RotationSpikesZPerc'};
writetable(MotTable,['.\HeadMotionPlots\' task{1} '_' sess(4:end) '\HeadMotionOutliers.xlsx']);
disp(MotTable);

end