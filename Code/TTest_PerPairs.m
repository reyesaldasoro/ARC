%% Clear all variables and close all figures
clear all
close all
clc
%% Read the files that have been stored in the current folder
if strcmp(filesep,'/')
    % Running in Mac    
    cd ('/Users/ccr22/Academic/GitHub/ARC/Code')
else
    % running in windows
    cd ('D:\Acad\GitHub\ARC\Code')
end
% Calculate the number of files, as of February 2019 there were 315 valid


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

%%
GroupsPresent   = unique (cumulativeStats(:,1));


%% Test per pairs

for caseBoxplot =4:9
    % Test all the metrics as described above
    ResultsTTest{caseBoxplot}(1:13,1:14)=nan;
    % Test all the pairs of groups
    for group1=GroupsPresent(1):GroupsPresent(end)-1
        for group2=group1+1:GroupsPresent(end)
            [h(group1,group2),ResultsTTest{caseBoxplot}(group1,group2),p]=ttest2(...
                cumulativeStats(cumulativeStats(:,1)==group1,caseBoxplot),...
                cumulativeStats(cumulativeStats(:,1)==group2,caseBoxplot));
        end
    end
end
