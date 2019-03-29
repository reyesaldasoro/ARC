%% Clear all variables and close all figures
clear all
close all
clc

%% Read the files that have been stored in the current folder
if strcmp(filesep,'/')
    % Running in Mac
    load('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Datasets/DataARC_Datasets_2019_02_27.mat')
    cd ('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Results')
    baseDir                             = 'Metrics_2019_03_14/metrics/';
else
    % running in windows
    load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2019_02_27.mat')
    cd ('D:\OneDrive - City, University of London\Acad\ARC_Grant\Results')
    baseDir                             = 'Metrics_2019_03_14/metrics/';
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
cumulativeStats         =[[] []];
cumulativeStatsTort     =[[] []];


for currentCase=1: numMetrics
    clear cent* min* max* rr cc relP* jet*
    % Load the current file to be displayed
    load(strcat(baseDir,MetricsDir(currentCase).name));
    qqq{currentCase}= MetricsDir(currentCase).name;
    %    cumulativeStats = [ cumulativeStats; [cell_metrics.Dist_um_s]' [nuclei_metrics.PositionR]' ];
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
    %10 nuclei_metrics.Dist_um_s
    %11 cell_metrics.Area_um_2
    %12 nuclei_metrics.Area_um_2
    %13 12 ./ 11 ratio of areas
    
    
    
    clear cumulativeStatsCurr
    % Filter the outliers
    try
        cumulativeStatsCurr(:,1:2)     =   ones(numCurrFrames,1)*[GroupCell(currentCase) currentCase];
        cumulativeStatsCurr(:,3)       =   [1:numCurrFrames]'  ;
        cumulativeStatsCurr(:,4)       =   imfilter([cell_metrics.Dist_um_s]',     [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,5)       =   imfilter([nuclei_metrics.PositionR]',   [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,6)       =   imfilter([nuclei_metrics.Min_MajAxis],  [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,7)       =   imfilter([nuclei_metrics.forkness]',    [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,8)       =   imfilter([cell_metrics.forkness]',      [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,9)       =   imfilter([cell_metrics.skelAlignment]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,10)      =   imfilter([nuclei_metrics.Dist_um_s]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,11)      =   imfilter([cell_metrics.Area_um_2]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,12)      =   imfilter([nuclei_metrics.Area_um_2]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,13)      =   cumulativeStatsCurr(:,12)./cumulativeStatsCurr(:,11);
        
    catch
        index1 =find(([nuclei_metrics.PositionR]));
        index2 =find(~isnan([nuclei_metrics.PositionR]));
        index3 =intersect(index1,index2);
        cumulativeStatsCurr(index3,4)       =   imfilter([cell_metrics.Dist_um_s]',     [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(:,5)            =   imfilter([nuclei_metrics.PositionR]',   [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,6)       =   imfilter([nuclei_metrics.Min_MajAxis],  [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,7)       =   imfilter([nuclei_metrics.forkness]',    [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,8)       =   imfilter([cell_metrics.forkness]',      [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,9)       =   imfilter([cell_metrics.skelAlignment]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,10)      =   imfilter([nuclei_metrics.Dist_um_s]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,11)      =   imfilter([cell_metrics.Area_um_2]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,12)      =   imfilter([nuclei_metrics.Area_um_2]', [0.25 0.5 0.25]','replicate');
        cumulativeStatsCurr(index3,13)      =   cumulativeStatsCurr(:,12)./cumulativeStatsCurr(:,11);
        
        
        
    end
    % cumulativeStatsCurr(:,5)       = imfilter(cumulativeStatsCurr(:,1),[0.25 0.5 0.25]','replicate');
    cumulativeStats                 = [ cumulativeStats; cumulativeStatsCurr];
    cumulativeStatsTort             = [ cumulativeStatsTort; [cumulativeStatsCurr(1,1:2) cell_metrics.pathTortuosity]]; 
end

labels={'group','case','time','Vel [um/s] (C)','Rel Position','Min/Maj','Forkness (N)','Forkness (C)','Skel Alignment','Vel [um/s] (N)','Area [um2] (C)','Area [um2] (N)','Area N/Area C'};
%% Results 2019_02_06
% As of 2019_02_06, this produces 10,232 rows, one for each instance of a
% cell with all its stats

%% remove the cases with NaNs (first of every track, no velocity is calculated)
cumulativeStats(isnan(cumulativeStats(:,4)),:)=[];

% As of 2019_02_06 there were 750 cases removed leaving 9483
%% remove the cases with speed above 0.49
cumulativeStats(cumulativeStats(:,4)>0.45,:)=[];

% As of 2019_02_06 there were 14 cases removed leaving 9468
%% remove the cases with skeletal alignment above 300
cumulativeStats(cumulativeStats(:,9)>200,:)=[];

% As of 2019_02_06 there were 197 cases removed leaving 9271
%% remove discarded cases
cumulativeStats(cumulativeStats(:,1)==0,:)=[];

% No cases discarded, all were removed previously
%% Display done in other files, just save in the GitHub folder

if strcmp(filesep,'/')
    % Running in Mac
    cd ('/Users/ccr22/Academic/GitHub/ARC/Code')
else
    % running in windows
    cd('D:\Acad\GitHub\ARC\Code')
end




filename =strcat('cumulativeStats_',datestr(date,'yyyy_mm_dd'));
save(filename,'cumulativeStats','labels','cumulativeStatsTort')

