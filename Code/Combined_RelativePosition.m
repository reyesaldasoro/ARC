%% Clear all variables and close all figures
clear all
close all
clc

%% Read the files that have been stored in the current folder
if strcmp(filesep,'/')
    % Running in Mac
    load('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Datasets/DataARC_Datasets_2019_05_03.mat')
    cd ('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Results')
    baseDir                             = 'Metrics_2019_04_25/metrics/';
else
    % running in windows
    load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2019_05_03.mat')
    cd ('D:\OneDrive - City, University of London\Acad\ARC_Grant\Results')
    baseDir                             = 'Metrics_2019_04_25/metrics/';
end
%% Calculate the number of files, as of February 2019 there were 315 valid
% files
MetricsDir                          = dir(strcat(baseDir,'*.mat'));
numMetrics                          = size(MetricsDir,1);


%%
figure
clf
hold on
q=nan(numMetrics,300);
for k=1: numMetrics
    % Load the current file to be displayed
    caseName                = MetricsDir(k).name;     
    currTrack               = strcat(baseDir,caseName); 
    load(currTrack);
    % Discard all cases with NaN, these are mainly when the cell is out of the field
    % of view, also when the centroid is a NaNs, it is recorded as 1 value, when it
    % should be 2
    nanVect                 = ~isnan([cell_metrics.Area]);
    timePoints              = numel(nanVect);
    % create a time and position vector 
    posVector               = k*ones(size(nanVect));
    timeVector              = 1:timePoints;
    plot3(timeVector,posVector,[nuclei_metrics.PositionR])
    q(k,1:timePoints)       = [nuclei_metrics.PositionR];
end
view(10,30)
grid on
axis tight

%% Calculate mean values, it has nans so mean will not work directly

q3                          = q;
q3(:,1)                     = q3(:,2);
%meanRelPos                  = sum(q3,1)./(sum(q3~=0));
%stdRelPos                   = sum(q3,1)./(sum(q3~=0));
for k = 1:150
    meanRelPos(k)       = mean(q3(~isnan(q3(:,k)),k));
    stdRelPos(k)        = std(q3(~isnan(q3(:,k)),k));  
end
%
figure
hold off
plot(q','color',0.8*[1 1 1])
hold on
plot(meanRelPos,'color',[1 0 0],'linewidth',3)
plot(meanRelPos+stdRelPos,'color',[1 0 0],'linewidth',1,'linestyle','--')
plot(meanRelPos-stdRelPos,'color',[1 0 0],'linewidth',1,'linestyle','--')
plot([1 150],[0 0],'k-')
grid on
xlabel('Time Frame','fontsize',20)
ylabel('Relative position of Nuclei','fontsize',20)
set(gcf,'Position',[33   412   967   386]);
axis([0 100 -0.8 0.8])
%% Calculate a vector for injured non-injured
%load('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Datasets/DataARC_Datasets_2018_03_20.mat')

%%
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


%%
%meanRelPos                  = sum(q3,1)./(sum(q3~=0));
%stdRelPos                   = sum(q3,1)./(sum(q3~=0));
q_Injured                    = q3(IsInjured==2,:);
q_NonInjured                 = q3(IsInjured==1,:);

for k = 1:150
    % Iterate over timepoints
        meanRelPos_NI(k)       = mean(q_NonInjured(~isnan(q_NonInjured(:,k)),k));
        stdRelPos_NI(k)        = std(q_NonInjured(~isnan(q_NonInjured(:,k)),k));        
        meanRelPos_I(k)       = mean(q_Injured(~isnan(q_Injured(:,k)),k));
        stdRelPos_I(k)        = std(q_Injured(~isnan(q_Injured(:,k)),k));
end
%%
figure
hold off
plot(meanRelPos_NI,'color',[1 0 0],'linewidth',3)
hold on
plot(meanRelPos_I,'color',[0 0 1],'linewidth',3)
plot(q','color',0.8*[1 1 1])
plot(meanRelPos_NI,'color',[1 0 0],'linewidth',3)
plot(meanRelPos_I,'color',[0 0 1],'linewidth',3)

plot(meanRelPos_NI+stdRelPos_NI,'color',[1 0 0],'linewidth',1,'linestyle','--')
plot(meanRelPos_NI-stdRelPos_NI,'color',[1 0 0],'linewidth',1,'linestyle','--')
plot(meanRelPos_I+stdRelPos_I,'color',[0 0 1],'linewidth',1,'linestyle','--')
plot(meanRelPos_I-stdRelPos_I,'color',[0 0 1],'linewidth',1,'linestyle','--')


hLegend = legend({'Non-injured','Injured'});


plot([1 150],[0 0],'k-')
grid on
xlabel('Time Frame','fontsize',20)
ylabel('Relative position of Nuclei','fontsize',20)
set(gcf,'Position',[33   412   967   386]);
axis([0 100 -0.8 0.8])

hLegend.String ={'Non-injured','Injured'};

