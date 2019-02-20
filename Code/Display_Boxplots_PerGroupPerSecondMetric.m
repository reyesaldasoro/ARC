%% Clear all variables and close all figures
clear all
close all

%% Read the current cumulative stats
% Find the files in the folder
CumulativeStats_Dir                     = dir('CumulativeStats_2019*');
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

%% Display ONE METRIC against the groups 1,2,3,... 14 AND A second metric divided by steps in boxplots
step = 5;
invStep = 1/step;
currentMetric_x           = 4;
currentMetric_y           = 5;
boxplot(cumulativeStats(:,currentMetric_y) ,{cumulativeStats(:,1),invStep/2+(invStep)*round(step*cumulativeStats(:,currentMetric_x))},'whisker',2)
grid on
xlabel({'Group';labels{currentMetric_x}},'fontsize',15)
ylabel(labels{currentMetric_y},'fontsize',20)
