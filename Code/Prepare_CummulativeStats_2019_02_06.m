%% Clear all variables and close all figures
clear all
close all

%% Read the files that have been stored in the current folder
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
%% Results 2019_02_06
% As of 2019_02_06, this produces 10,232 rows, one for each instance of a
% cell with all its stats

%% remove the cases with NaNs (first of every track, no velocity is calculated)
cummulativeStats(isnan(cummulativeStats(:,4)),:)=[];

% As of 2019_02_06 there were 750 cases removed leaving 9483
%% remove the cases with speed above 0.49
cummulativeStats(cummulativeStats(:,4)>0.45,:)=[];

% As of 2019_02_06 there were 14 cases removed leaving 9468
%% remove the cases with skeletal alignment above 300
cummulativeStats(cummulativeStats(:,9)>200,:)=[];

% As of 2019_02_06 there were 197 cases removed leaving 9271
%% remove discarded cases
cummulativeStats(cummulativeStats(:,1)==0,:)=[];

% No cases discarded, all were removed previously
%% Display done in other files, just save

save(strcat('CummulativeStats_',datestr(date,'yyyy_mm_dd')),cummulativeStats)

