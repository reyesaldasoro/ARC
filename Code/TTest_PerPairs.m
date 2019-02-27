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
labels={'group','case','time','Dist [um/s]','Rel Position','Min/Maj','Forkness (N)','Forkness (C)','Skel Alignment'};

%%
GroupsPresent   = unique (cumulativeStats(:,1));


%% Test per pairs

for caseBoxplot =4:9
    % Test all the metrics as described above
    ResultsTTest{caseBoxplot}(1:13,1:14)= nan;
    h_TTest{caseBoxplot}(1:13,1:14)     = nan;
    
    % Test all the pairs of groups
    for group1=GroupsPresent(1):GroupsPresent(end)-1
        for group2=group1+1:GroupsPresent(end)
            [h_TTest{caseBoxplot}(group1,group2),ResultsTTest{caseBoxplot}(group1,group2),p]=ttest2(...
                cumulativeStats(cumulativeStats(:,1)==group1,caseBoxplot),...
                cumulativeStats(cumulativeStats(:,1)==group2,caseBoxplot));
        end
    end
end


%% ANOVAS
% Lamin mutants
% Wt (4) vs. A-/-(12), B1-/-(13), B2-/-(14)
% (12) vs. (13)
% (12) vs (14)
% (13) vs (14)
% Lbr mutant
% WT (4) vs. -/- lbr mutant (5)
% Simple T test (unpaired)
% 
% Wt (4) vs. lbr+/-(16) vs lbr-/-(5)
% Ordered comparison test (not ANOVA) as expect a gene dosage response where (16) is an intermediate between (4) and (5). I’m not sure what the correct test would be.
%  
% Uninjured vs. injured
%  
% Uninjured (1) vs. (2) injured
% T-test (unpaired)
%  
% (2) injured B2 reporter background vs (6) injured H2Bcerulean background
% As a control comparison for effect of different transgenes (B2 vs H2B)
% T-test (unpaired)
%  
% (2) injured B2 reporter vs. (7) injured B2 reporter DMSO 
% As a control to confirm that DMSO is not having an overt effect
% T-test (unpaired)

% ANOVA1 for a series of groups
currentMetric           = 7;
groupsToSelect          = [5 12 13 14];
selectionPoints         = ismember(cumulativeStats(:,1),groupsToSelect);
[p,t,stats]= anova1(cumulativeStats(selectionPoints,currentMetric),cumulativeStats(selectionPoints,1),'off');


% Perform a multiple comparison of the group means.
figure(currentMetric)
[c,m,h,nms] = multcompare(stats);
ylabel(labels{currentMetric})
%%
% Anova1 for a series of groups, some of these combined
% currentMetric           = 4;
% groups1                 = [6 ];
% groups2                 = [12 13 14];
% 
% selectionPoints1         = ismember(cumulativeStats(:,1),groups1);
% selectionPoints2         = ismember(cumulativeStats(:,1),groups2);
% group1                  = cumulativeStats(selectionPoints1,currentMetric);
% group2                  = cumulativeStats(selectionPoints2,currentMetric);
% 
% d= anova1([group1;group2],[ones(size(group1));2*ones(size(group2))]);



