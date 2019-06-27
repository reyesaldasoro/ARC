%% Clear all variables and close all figures
clear all
close all

%% Read the current cumulative stats
% Find the files in the folder
CumulativeStats_Dir                     = dir('CumulativeStats_2019*.mat');
% Take the latest one and load the stats
load(CumulativeStats_Dir(end).name,'cumulativeStats','labels')
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
%labels={'group','case','time','Dist [um/s]','Rel Position','Min/Maj','Forkness (N)','Forkness (C)','Skel Alignment'};

%% Display ONE METRIC against ALL the groups 1,2,3,... 14
currentMetric_x           = 1;
currentMetric_y           = 13;
%figure
boxplot(cumulativeStats(:,currentMetric_y),cumulativeStats(:,currentMetric_x),'whisker',1)
grid on
xlabel(labels{currentMetric_x},'fontsize',20)
ylabel(labels{currentMetric_y},'fontsize',20)
 %axis([0.5 5.5 -0.75 0.6])
 %axis([0.25 10.75 -0.75 0.65])

%% Display ONE METRIC against A SUBSET of the groups 1,2,3,... 14
currentMetric_x             =  1;
currentMetric_y             = 13;
subsetGroups                = [1 3 6 14 16];
%figure
indexGroups                 = ismember(cumulativeStats(:,currentMetric_x),subsetGroups);
boxplot(cumulativeStats(indexGroups,currentMetric_y),cumulativeStats(indexGroups,currentMetric_x),'whisker',1)
grid on
xlabel(labels{currentMetric_x},'fontsize',20)
ylabel(labels{currentMetric_y},'fontsize',20)
 %axis([0.5 5.5 -0.75 0.6])
 %axis([0.25 10.75 -0.75 0.65])
