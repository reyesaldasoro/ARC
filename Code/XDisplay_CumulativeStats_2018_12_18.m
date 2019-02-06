%% Clear all variables and close all figures
clear all
close all

%% Read the current cumulative stats
% Find the files in the folder
CumulativeStats_Dir                     = dir('CumulativeStats_2019*');
% Take the latest one and load the stats
load(CumulativeStats_Dir(end).name)
%%

    % calculate the stats of the current video In the Following Order:
    % 1 GROUP, 1-non injured, etc
    % 2 CASE, i.e. the order of the files
    % 3 Time point
    % 4 cell_metrics.Dist_um_s
    % 5 nuclei_metrics.PositionR
    % 6 nuclei_metrics.Min_MajAxis
    % 7 nuclei_metrics.forkness
    % 8 cell_metrics.forkness
    % 9 cell_metrics.skelAlignment
%% Prepare labels
    
labels={'group','case','time','Dist [um/s]','Rel Position','Min/Maj','Forkness (N)','Forkness (C)','Skel Alignment'};

%% Display ... rel position against the velocity divided by steps in boxplots
step = 10;
invStep = 1/step;
currentMetric_x           = 4;
currentMetric_y           = 5;
figure
boxplot(cummulativeStats(:,5),invStep/2+(invStep)*round(step*cummulativeStats(:,4)),'whisker',1)
grid on
xlabel(labels{currentMetric_x},'fontsize',20)
ylabel(labels{currentMetric_y},'fontsize',20)
 %axis([0.5 5.5 -0.75 0.6])
 %axis([0.25 10.75 -0.75 0.65])
%% rel position against the velocity divided by steps 
figure
plot(invStep/2+(invStep)*round(step*cummulativeStats(:,5)),cummulativeStats(:,2),'x')
grid on
xlabel('Inst. velocity','fontsize',20)
ylabel('Relative position of Nuclei','fontsize',20)

%% Display ... rel position against the velocity divided by steps AND injured/non-injured in boxplots
step = 10;
invStep = 1/step;

caseBoxplot = 5;
% figure
% boxplot(cummulativeStats(:,caseBoxplot) ,{invStep/2+(invStep)*round(step*cummulativeStats(:,4)),cummulativeStats(:,1)})
% grid on
% xlabel(labels{4},'fontsize',20)
% ylabel(labels{caseBoxplot},'fontsize',20)


figure
boxplot(cummulativeStats(:,caseBoxplot) ,{cummulativeStats(:,1),invStep/2+(invStep)*round(step*cummulativeStats(:,4))})
grid on
xlabel(labels{4},'fontsize',20)
ylabel(labels{caseBoxplot},'fontsize',20)



%%
caseBoxplot =9;
figure(caseBoxplot)
boxplot(cummulativeStats(:,caseBoxplot),cummulativeStats(:,1))
grid on
ylabel(labels{caseBoxplot})


%caseBoxplot =6;
group1 = 5;
group2 = 6;

for group1=1:5
    for group2=group1+1:6
        [h,t(group1,group2),p]=ttest2(cummulativeStats(cummulativeStats(:,1)==group1,caseBoxplot),cummulativeStats(cummulativeStats(:,1)==group2,caseBoxplot));
    end
end