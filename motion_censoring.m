% Motion censoring according to Power et al., 2012 / Siegel et al., 2014
%
% Framewise Displacement (FD)
% FD = |Delta dix| + |Delta diy| + |Delta diz| + |Delta alphai| + |Delta betai| + |Delta gammai|

% Christian Baeuchl | last edit 01.06.2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear variables, close windows and clear the command line
clc,clear,close all
tic

path             = 'path to your preprocessed data';
scriptpath       = 'path to your code files';
returnpath       = 'path to location where your outlier files will be stored (folder must be created in advance)';
skip_plotting    = 0;
additional_plots = 0;

list = sort(importdata('.\Participant_lists\List.txt'));

functfold   = 'folder with your functional images';
session     = 'BL';
FDcrit      = 0.9;
brainradius = 50; % mm
n = 1;

%%
fprintf(1,'\nCalculating Framewise Displacement (FD)\n')
fprintf(1,'_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n');
fprintf(1,'\n');
fprintf(1,'%d data sets in total\n',size(list,1))
fprintf(1,'Processing data set #: ' );

for i = 1:length(list)

    %% calculating FD
    
    if i>1
        for b=0:log10(i-1)
            fprintf(1,'\b');  % Deleting last displayed number
        end
    end
    fprintf(1,'%d',i);

    cd([path num2str(list(i)) filesep functfold ])
    current_path = [path num2str(list(i)) filesep functfold filesep];
    rpname = ls([current_path 'Vxrp*.txt']);
    if isempty(rpname)
        warning('ID: %s is missing the Vxrp file. Skipping subject\n ',num2str(list(i)))
        continue
    end
    rp      = load([current_path rpname]);
    rpname2 = ls([current_path 'rp*.txt']);
    rp2     = load([current_path rpname2]);

    FD = abs(rp(:,7)) + abs(rp(:,8)) + abs(rp(:,9)) + abs(rp(:,10))*brainradius + abs(rp(:,11))*brainradius + abs(rp(:,12))*brainradius;
    FDMEAN = mean(FD);
    FDMEAN = round(FDMEAN,3);
    FDMEDIAN = median(FD);
    FDMEDIAN = round(FDMEDIAN,3);
    RMS = rms([(rp(:,1));(rp(:,2));(rp(:,3));(rp(:,4)*brainradius);(rp(:,5)*brainradius);(rp(:,6)*brainradius)]);
    RMS = round(RMS,3);
    FDRMS = rms([(rp(:,7));(rp(:,8));(rp(:,9));(rp(:,10)*brainradius);(rp(:,11)*brainradius);(rp(:,12)*brainradius)]);
    FDRMS = round(FDRMS,3);
    tempplot = zeros(length(FD),1);

    %%
    if any(FD(FD>FDcrit))
        tempplot(FD>FDcrit) = 1;
        indx = find(FD>FDcrit);
        tempmask = zeros(length(FD),numel(indx));
        for j = 1:numel(indx)
            tempmask(indx(j),j) = 1;
        end
        checksum = sum(tempplot == sum(tempmask,2));
        if checksum ~= length(FD)
            warning('Check FD Regressors of subject/session: %s\n',num2str(list(i)))
        end
        rpnew = [tempmask,rp];
        save([returnpath filesep 'FD_Vxrp_' num2str(list(i)) '.txt'],'rpnew','-ascii','-double');
        save([current_path filesep 'FD_Vxrp_' num2str(list(i)) '.txt'],'rpnew','-ascii','-double');
        rpname2 = ls([current_path 'rp_*.txt']);
        rp2 = load([current_path rpname2]);
        rpnew2 = [tempmask,rp2];
        save([returnpath filesep   'FD_rp_' num2str(list(i)) '.txt'],'rpnew2','-ascii','-double');
        save([current_path filesep 'FD_rp_' num2str(list(i)) '.txt'],'rpnew2','-ascii','-double');
    else
        save([returnpath filesep   'FD_Vxrp_' num2str(list(i)) '.txt'],'rp','-ascii','-double');
        save([current_path filesep 'FD_Vxrp_' num2str(list(i)) '.txt'],'rp','-ascii','-double');
        save([returnpath filesep   'FD_rp',num2str(list(i)) '.txt'],'rp2','-ascii','-double');
        save([current_path filesep 'FD_rp',num2str(list(i)) '.txt'],'rp2','-ascii','-double');
    end

    %% compile list
    ID{n} = num2str(list(i));
    meanFD(n) = FDMEAN;
    medianFD(n) = FDMEDIAN;
    if ~exist('indx','var')
        numcritFD(n) = 0;
        indx = [];
        percFD(n) = 0;
    else
        numcritFD(n) = numel(indx);
        percFD(n)    = (100*numel(indx))/length(FD);
    end
    n = n + 1;

    %% plotting
    if ~skip_plotting
        map = [0.25 0.65 1]; % color map
        figure;
        set(gcf,'color','w')
        subplot(5,1,1:4)
        ar = area(FD);
        ar.FaceColor = map;
        colormap(map)
        hold on
        title({['\fontsize{15}' num2str(list(i)) '_' session], ['\fontsize{13}RMS movement: ' num2str(RMS) ' mm ' '| FD RMS movement: ' num2str(FDRMS) ' mm ' '| FD(mean): ' num2str(FDMEAN) ' mm ' '| # Scans > FD crit.: ' num2str(numel(indx))]});
        ylabel('FD (mm)')
        plot(FDcrit*ones(1,length(FD)),'--','color','k'); % line for outliers
        dim = [0.35 0.1 .85 .75];
        axis([1 length(FD) 0 max([max(FD),FDcrit])+0.2])
        set(gca,'XTickLabel',[]);
        subplot(5,1,5)
        plot(tempplot,'-','color',map)
        axis([1 length(FD) 0 1])
        set(gca,'YTickLabel',[]);
        xlabel('scans')
        set(gcf,'PaperOrientation','landscape');
        hold off
        print([returnpath filesep num2str(list(i)) '_FD'],'-fillpage','-dpdf');
        close(gcf)
    end

    %% clean up
    clear current_path indx indx2 rp rp2 rpname rpname2 rpnew rpnew2 FD FDRMS FDMEAN FDMEDIAN tempmask tempplot

end

if (n-1) ~= i
    warning('Not all subjects from list were processed')
end

percFD = round(percFD,3);
%% create and print table
T = table(ID',meanFD',medianFD',numcritFD',percFD');
T.Properties.VariableNames = {'ID','FD_mean','FD_median','FD_crit','FD_crit_perc'};
writetable(T,fullfile(returnpath,['FD_stats_' functfold(1:end-4) '_' session '_' datestr(now,'yymmdd'),'.xlsx']),'Sheet','all');
disp(T); %,T2

if additional_plots
    %%
    filename = cfg_getfile('FPListRec',fullfile(returnpath),'^FD_stats_FOO_231013.xlsx');
    data = readtable(filename{1},'Sheet','all');
    AVG      = 'FD_mean'; % 'FD_median';
    data_PLA = data(strcmp(data.('Intervention'),'PLACEBO'),AVG);
    data_LDO = data(strcmp(data.('Intervention'),'LDOPA'),AVG);
    data_PLS = data(strcmp(data.('Intervention_Order'),'PLS'),AVG);
    data_LDS = data(strcmp(data.('Intervention_Order'),'LDS'),AVG);
    data_YA  = data(strcmp(data.('AgeGroup'),'YA'),AVG);
    data_OA  = data(strcmp(data.('AgeGroup'),'OA'),AVG);
    data_S1  = data(data.('MRI_Session') == 1,AVG);
    data_S2  = data(data.('MRI_Session') == 2,AVG);
    plotdata = data.('FD_mean'); 
    plotdataGRP1 = data_PLA{:,1}; %data_PLS{:,1}; %data_YA{:,1}; % data_S1{:,1}; % data_PLA{:,1}; 
    plotdataGRP2 = data_LDO{:,1}; %data_LDS{:,1}; %data_OA{:,1}; % data_LDO{:,1}; % data_S2{:,1};
    plotdataname = 'All data';
    plotdataGRP1name = 'PLACEBO';%'PLS'; %'YA'; %'Session 1'; %'PLACEBO'; 
    plotdataGRP2name = 'L-DOPA';%'LDS'; %'OA'; %'Session 2'; %'L-DOPA'; 
    plotAVGname = 'FD (mean)'; % 'FD (median)';

    %
    figure
    h = boxplot(plotdata,'Notch','on','Labels',{'All subs (S1 & S2)'},'Whisker',1);
    set(h, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h,'Color','k')
    ylabel('FD(mean)')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c = [0.3 0.3 0.3];
    x=ones(length(plotdata)).*(1+(rand(length(plotdata))-0.5)/15);
    f1=scatter(x(:,1),plotdata(:,1),[],c,'filled');f1.MarkerFaceAlpha = 0.4;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:round(max(plotdata),1)])
    daspect([1 0.75 1])
    ylim([0,round(max(plotdata),1))
    hold off

    %
    figure
    h2 = boxplot(plotdataGRP1,'Notch','on','Labels',{'Session 1'},'Whisker',1);
    set(h2, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h2,'Color','k')
    ylabel('FD(mean)')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c =  [1 0.5 0.5]; %[0.9290, 0.6940, 0.1250]; %[0.8500 0.3250 0.0980]; % % [0.6350 0.0780 0.1840];
    x2=ones(length(plotdataGRP1)).*(1+(rand(length(plotdataGRP1))-0.5)/15);
    f2=scatter(x2(:,1),plotdataGRP1(:,1),[],c,'filled');f2.MarkerFaceAlpha = 0.6;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:round(max(plotdataGRP1),1)])
    ylim([0,round(max(plotdataGRP1),1)])
    daspect([1 0.75 1])
    hold off

    %
    figure
    h3 = boxplot(plotdataGRP2,'Notch','on','Labels',{'Session 2'},'Whisker',1);
    set(h3, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h3,'Color','k')
    ylabel('FD(mean)')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c =  [0.25 0.5 1]; %[0 0.4470 0.7410]; %[0.3010 0.7450 0.9330]; %
    x3=ones(length(plotdataGRP2)).*(1+(rand(length(plotdataGRP2))-0.5)/15);
    f3=scatter(x3(:,1),plotdataGRP2(:,1),[],c,'filled');f3.MarkerFaceAlpha = 0.6;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:round(max(plotdataGRP1),1)])
    ylim([0,round(max(plotdataGRP1),1)])
    daspect([1 0.75 1])
    hold off

    %%
    figure
    subplot(1,3,1)
    h = boxplot(plotdata,'Notch','on','Labels',{plotdataname},'Whisker',1);
    set(h, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h,'Color','k')
    ylabel('FD(mean)')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c = [0.3 0.3 0.3];
    x = ones(length(plotdata)).*(1+(rand(length(plotdata))-0.5)/15);
    f1 = scatter(x(:,1),plotdata(:,1),[],c,'filled');f1.MarkerFaceAlpha = 0.4;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:round(max(plotdata),1)])
    daspect([1 0.75 1])
    ylim([0,round(max(plotdata),1)])
    hold off

    subplot(1,3,2)
    h2 = boxplot(plotdataGRP1,'Notch','on','Labels',{plotdataGRP1name},'Whisker',1);
    set(h2, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h2,'Color','k')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c =  [1 0.5 0.5]; %[0.9290, 0.6940, 0.1250]; %[0.8500 0.3250 0.0980]; % % [0.6350 0.0780 0.1840];
    x2=ones(length(plotdataGRP1)).*(1+(rand(length(plotdataGRP1))-0.5)/15);
    f2=scatter(x2(:,1),plotdataGRP1(:,1),[],c,'filled');f2.MarkerFaceAlpha = 0.6;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:round(max(plotdata),1)])
    ylim([0,round(max(plotdata),1)])
    daspect([1 0.75 1])
    hold off

    subplot(1,3,3)
    h3 = boxplot(plotdataGRP2,'Notch','on','Labels',{plotdataGRP2name},'Whisker',1);
    set(h3, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h3,'Color','k')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c =  [0.25 0.5 1]; %[0 0.4470 0.7410]; %[0.3010 0.7450 0.9330]; %
    x3=ones(length(plotdataGRP2)).*(1+(rand(length(plotdataGRP2))-0.5)/15);
    f3=scatter(x3(:,1),plotdataGRP2(:,1),[],c,'filled');f3.MarkerFaceAlpha = 0.6;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:round(max(plotdata),1)])
    ylim([0,round(max(plotdata),1)])
    daspect([1 0.75 1])
            set(gcf,'PaperOrientation','landscape');
    hold off



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if additional_plots
    [~,~,fulldata] = xlsread(char(filename),'Sheet3');
    data = xlsread(char(filename),'Sheet3');
    plotdata = data(:,1);
    plotdata_A = data(1:77,1);
    plotdata_B = data(78:end,1);

    %
    figure
    h4 = boxplot(plotdata_A,'Notch','on','Labels',{'Placebo'},'Whisker',1);
    set(h4, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h4,'Color','k')
    ylabel('FD(mean)')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c = [0 0.6 0];
    x4=ones(length(plotdata_A)).*(1+(rand(length(plotdata_A))-0.5)/15);
    f4=scatter(x4(:,1),plotdata_A(:,1),[],c,'filled');f4.MarkerFaceAlpha = 0.6;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:1.3])
    ylim([0,1.2])
    daspect([1 0.75 1])
    hold off

    %
    figure
    h5 = boxplot(plotdata_B,'Notch','on','Labels',{'L Dopa'},'Whisker',1);
    set(h5, 'DefaultTextFontSize', 30)
    set(findobj(gca,'type','line'),'linew',1.5)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(h5,'Color','k')
    ylabel('FD(mean)')
    set(lines, 'Color', 'r')
    set(gcf,'color','w')
    hold on
    c = [0.7 0.4 1];
    x5=ones(length(plotdata_B)).*(1+(rand(length(plotdata_B))-0.5)/15);
    f5=scatter(x5(:,1),plotdata_B(:,1),[],c,'filled');f5.MarkerFaceAlpha = 0.6;
    set(gca, 'ActivePositionProperty', 'position', 'FontWeight','bold','FontSize',18)
    set(gca,'ytick',[0:0.1:1.3])
    ylim([0,1.2])
    daspect([1 0.75 1])
    hold off

    %% simple group comparison

    [hi,p,ci,stats] = ttest(plotdataGRP1,plotdataGRP2,'Tail','both');
    [p,hi,stats] = signrank(plotdata_A,plotdata_B, 'tail','both');
end
toc
