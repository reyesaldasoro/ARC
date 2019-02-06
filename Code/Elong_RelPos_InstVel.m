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
%% Calculate a vector for injured non-injured
%load('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Datasets/DataARC_Datasets_2018_04_17.mat')
%load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2018_04_17.mat')

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
    IsInjured(k) = DataARC_Datasets{caseNumber,21};
end



%% read all files and accumulate the stats
cumulativeStats=[[] []];
for k=1: numMetrics
    clear cent* min* max* rr cc relP* jet*
    clear cumulativeStatsCurr
    % Load the current file to be displayed
    load(strcat(baseDir,MetricsDir(k).name));
    %    cumulativeStats = [ cumulativeStats; [cell_metrics.Dist_um_s]' [nuclei_metrics.PositionR]' ];
    numCurrFrames               = size(nuclei_metrics,2);
    
    
    
    % calculate the stats of the current video
    %    cumulativeStatsCurr        = [[cell_metrics.Dist_um_s]'  max([nuclei_metrics.Min_MajAxis])-[nuclei_metrics.Min_MajAxis]' ones(numCurrFrames,1)*[IsInjured(k) k] [1:numCurrFrames]'  ];
    
    try
        % If the cell does not start in time 1, there may be errors in the
        % allocation of values to cumulativeStatsCurr
        cumulativeStatsCurr        = [[cell_metrics.Dist_um_s]'  [nuclei_metrics.Min_MajAxis]' [nuclei_metrics.PositionR]' ones(numCurrFrames,1)*[IsInjured(k) k] [1:numCurrFrames]'  ];
    catch
        index1 =find(([nuclei_metrics.PositionR]));
        index2 =find(~isnan([nuclei_metrics.PositionR]));
        index3 =intersect(index1,index2);
        cumulativeStatsCurr(index3,1)  =    [cell_metrics.Dist_um_s]' ;
        cumulativeStatsCurr(index3,2)  =    [nuclei_metrics.Min_MajAxis]' ;
        cumulativeStatsCurr(:,3)  =    [nuclei_metrics.PositionR]';
        cumulativeStatsCurr(:,4:5)=    ones(numCurrFrames,1)*[IsInjured(k) k];
        cumulativeStatsCurr(:,6)  =    [1:numCurrFrames]' ;
        


    
    
    end
    
    cumulativeStatsCurr(1,1:2) = cumulativeStatsCurr(2,1:2);
    % Low pass filter results
    cumulativeStatsCurr(:,6)   = imfilter(cumulativeStatsCurr(:,1),[0.25 0.5 0.25]','replicate'); % Speed
    cumulativeStatsCurr(:,7)   = imfilter(cumulativeStatsCurr(:,2),[0.25 0.5 0.25]','replicate'); % Elongation
    cumulativeStatsCurr(:,8)   = imfilter(cumulativeStatsCurr(:,3),[0.25 0.5 0.25]','replicate'); % Position
   
%    cumulativeStats            = [ cumulativeStats; [cell_metrics.Dist_um_s]' [nuclei_metrics.PositionR]' ones(numCurrFrames,1)*[IsInjured(k) k] [1:numCurrFrames]'  ];
    cumulativeStats            = [ cumulativeStats; cumulativeStatsCurr];
%        [cell_metrics.Dist_um_s]' [nuclei_metrics.PositionR]' ones(numCurrFrames,1)*[IsInjured(k) k] [1:numCurrFrames]'  ];

end

%% remove the cases with NaNs (first of every track, no velocity is calculated)
cumulativeStats(isnan(cumulativeStats(:,1)),:)=[];
cumulativeStats(isnan(cumulativeStats(:,3)),:)=[];
%% remove discarded cases
cumulativeStats(cumulativeStats(:,3)==0,:)=[];

%% remove the cases with speed above 0.49
cumulativeStats(cumulativeStats(:,1)>0.44,:)=[];


%%
cumulativeStats(cumulativeStats(:,6)>0.45,:)=[];



%%

casePlots=[6 7 8;6 8 7;7 6 8; 7 8 6;8 6 7; 8 7 6 ];
selectCase =2;

step1       =10;
invStep1    = 1/step1;


select1                      = invStep1/2+(invStep1)*round(step1*cumulativeStats(:,casePlots(selectCase,1)));
select1u                     = unique(select1);
select1u (isnan(select1u))    = [];
select_1w                    = min(diff(select1u));
numItems_1                   = numel(select1u);

step2       =10;
invStep2    = 1/step2;


select2                      = invStep2/2+(invStep2)*round(step2*cumulativeStats(:,casePlots(selectCase,2)));
select2u                     = unique(select2);
select2u (isnan(select2u))    = [];
select_2w                    = min(diff(select2u));
numItems_2                   = numel(select2u);


figure

for k1=1:numItems_1
    select1_Curr                = select1==select1u(k1);
    for k2= 1:numItems_2
        select2_Curr                = select2==select2u(k2);
        
        y = quantile(cumulativeStats(select1_Curr&select2_Curr,casePlots(selectCase,3)),[0.02 0.25  0.50  0.75 0.975]);
        boxPlot3D_old(select1u(k1),select2u(k2),y,select_1w,select_2w);
    end
end
rotate3d on
view(3)
grid on
axis tight
% casePlots=[6 7 8;6 8 7;7 6 8; 7 8 6;8 6 7; 8 7 6 ];

switch selectCase
    case 1
        xlabel('Inst. velocity','fontsize',20)
        ylabel('Min/Maj Axis Nuclei','fontsize',20)
        zlabel('Relative position of Nuclei','fontsize',20)
    case 2
        xlabel('Inst. velocity','fontsize',20)
        zlabel('Min/Maj Axis Nuclei','fontsize',20)
        ylabel('Relative position of Nuclei','fontsize',20)
    case 3
        ylabel('Inst. velocity','fontsize',20)
        xlabel('Min/Maj Axis Nuclei','fontsize',20)
        zlabel('Relative position of Nuclei','fontsize',20)
    case 4
        zlabel('Inst. velocity','fontsize',20)
        xlabel('Min/Maj Axis Nuclei','fontsize',20)
        ylabel('Relative position of Nuclei','fontsize',20)
    case 5
        ylabel('Inst. velocity','fontsize',20)
        zlabel('Min/Maj Axis Nuclei','fontsize',20)
        xlabel('Relative position of Nuclei','fontsize',20)
    case 6
        zlabel('Inst. velocity','fontsize',20)
        ylabel('Min/Maj Axis Nuclei','fontsize',20)
        xlabel('Relative position of Nuclei','fontsize',20)
end