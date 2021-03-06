%% Clear all variables and close all figures
clear all
close all




%% Read the current cumulative stats
% Find the files in the folder
CumulativeStats_Dir                     = dir('CumulativeStats_2019*.mat');
% Take the latest one and load the stats
load(CumulativeStats_Dir(end).name,'cumulativeStats','labels','cumulativeStatsTort')
%% Current order of the stats:
    % 1 GROUP, 1-non injured, etc
    % 2 CASE, i.e. the order of the files
    % 3 Time point
    % 4 cell_metrics.Dist_um_s
    % 5 nuclei_metrics.PositionR
    % 6 nuclei_metrics.Min_MajAxis
    % 7 nuclei_metrics.forkness
    % 8 cell_metrics.forkness
    % 9 cell_metrics.skelAlignment
% And Labels:
labels={'group','case','Path Tortuosity','Dist [um/s]','Rel Position','Min/Maj','Forkness (N)','Forkness (C)','Skel Alignment'};

%% Display ONE METRIC against the groups 1,2,3,... 14
currentMetric_x           = 1;
currentMetric_y           = 3;
%figure
%boxplot(cumulativeStats(:,currentMetric_y),cumulativeStats(:,currentMetric_x),'whisker',1)
boxplot(cumulativeStatsTort(:,currentMetric_y),cumulativeStatsTort(:,currentMetric_x),'whisker',1)

grid on
xlabel(labels{currentMetric_x},'fontsize',20)
ylabel(labels{currentMetric_y},'fontsize',20)
 %axis([0.5 5.5 -0.75 0.6])
 %axis([0.25 10.75 -0.75 0.65])
set(gca','yscale','log')