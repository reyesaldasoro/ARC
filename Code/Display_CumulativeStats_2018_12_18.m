%% Clear all variables and close all figures
clear all
close all

%% Read the current cumulative stats
% Find the files in the folder
CumulativeStats_Dir                     = dir('CumulativeStats_2019*');
% Take the latest one and load the stats
load(CumulativeStats_Dir(end).name)
%%
if strcmp(filesep,'/')
    % Running in Mac
    load('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Datasets/DataARC_Datasets_2019_02_01.mat')
    cd ('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Results')
    baseDir                             = 'Metrics_2019_02_01/metrics/';
else
    % running in windows
    load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2019_02_01.mat')
    cd ('D:\OneDrive - City, University of London\Acad\ARC_Grant\Results')
    baseDir                             = 'Metrics_2019_02_01/metrics/';
end
% Calculate the number of files, as of February 2019 there were 315 valid
% files
MetricsDir                          = dir(strcat(baseDir,'*.mat'));
numMetrics                          = size(MetricsDir,1);

%% Calculate a vector for injured non-injured and groups
% This will be used to create groups for 
GroupCell(numMetrics,1)=0;
%This loop is necessary as the order of the files and the DataARC may not
%be the same
for k=1: numMetrics
    % Load the current file to be displayed
    caseName                = MetricsDir(k).name;
    currTrack               = strcat(baseDir,caseName);
    caseNumber              = 0;
    for counterCase         = 1:size(DataARC_Datasets,1)
        if strcmp(DataARC_Datasets{counterCase,2}(12:end-1),caseName(1:end-19))
            caseNumber      = counterCase;
        end
    end
    GroupCell(k,1) = DataARC_Datasets{caseNumber,21};
end



%% read all files and accummulate the stats
cummulativeStats=[[] []];

for currentCase=1: numMetrics
    clear cent* min* max* rr cc relP* jet*
    % Load the current file to be displayed
    load(strcat(baseDir,MetricsDir(currentCase).name));
    %    cummulativeStats = [ cummulativeStats; [cell_metrics.Dist_um_s]' [nuclei_metrics.PositionR]' ];
    numCurrFrames               = size(nuclei_metrics,2);
   
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
    clear cummulativeStatsCurr
    % Filter the outliers
    try
        cummulativeStatsCurr(:,1:2)     =   ones(numCurrFrames,1)*[GroupCell(currentCase) currentCase];
        cummulativeStatsCurr(:,3)       =   [1:numCurrFrames]'  ;
        cummulativeStatsCurr(:,4)       =   imfilter([cell_metrics.Dist_um_s]',     [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(:,5)       =   imfilter([nuclei_metrics.PositionR]',   [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(:,6)       =   imfilter([nuclei_metrics.Min_MajAxis],  [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(:,7)       =   imfilter([nuclei_metrics.forkness]',    [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(:,8)       =   imfilter([cell_metrics.forkness]',      [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(:,9)       =   imfilter([cell_metrics.skelAlignment]', [0.25 0.5 0.25]','replicate');
    catch
        index1 =find(([nuclei_metrics.PositionR]));
        index2 =find(~isnan([nuclei_metrics.PositionR]));
        index3 =intersect(index1,index2);
        cummulativeStatsCurr(index3,4)       =   imfilter([cell_metrics.Dist_um_s]',     [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(:,5)       =   imfilter([nuclei_metrics.PositionR]',   [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(index3,6)       =   imfilter([nuclei_metrics.Min_MajAxis],  [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(index3,7)       =   imfilter([nuclei_metrics.forkness]',    [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(index3,8)       =   imfilter([cell_metrics.forkness]',      [0.25 0.5 0.25]','replicate');
        cummulativeStatsCurr(index3,9)       =   imfilter([cell_metrics.skelAlignment]', [0.25 0.5 0.25]','replicate');
    end
    % cummulativeStatsCurr(:,5)       = imfilter(cummulativeStatsCurr(:,1),[0.25 0.5 0.25]','replicate');
    cummulativeStats                = [ cummulativeStats; cummulativeStatsCurr];

end

labels={'group','case','time','Dist [um/s]','Rel Position','Min/Maj','Forkness (N)','Forkness (C)','Skel Alignment'}
%% remove the cases with NaNs (first of every track, no velocity is calculated)
cummulativeStats(isnan(cummulativeStats(:,4)),:)=[];

%% remove the cases with speed above 0.49
cummulativeStats(cummulativeStats(:,4)>0.45,:)=[];

%% remove the cases with skeletal alignment above 300
cummulativeStats(cummulativeStats(:,9)>200,:)=[];


%% remove discarded cases
cummulativeStats(cummulativeStats(:,1)==0,:)=[];

%% Display ... rel position against the velocity divided by steps in boxplots
step = 10;
invStep = 1/step;

figure
boxplot(cummulativeStats(:,5),invStep/2+(invStep)*round(step*cummulativeStats(:,4)),'whisker',1)
grid on
xlabel('Inst. velocity','fontsize',20)
ylabel('Relative position of Nuclei','fontsize',20)
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