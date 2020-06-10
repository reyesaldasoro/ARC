clear all; close all;clc
%%
if strcmp(filesep,'/')
    % Running in Mac
    cd '/Users/ccr22/Academic/ARC_Grant/ARC'   % Run in PowerMac
    %cd /Volumes/TOSHIBA/ARC/                    % Run in Laptop
    load('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Datasets/DataARC_Datasets_2019_05_03.mat')
else
    % running in windows
    cd ('D:\Acad\ARC\ARC')
%    load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2019_05_03.mat')
    load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2020_03_23.mat')

end
% Read all folders with data
dir0            = dir ('LSM*');
numFolders_0    = size(dir0,1);
%%
baseFolder_num      = 55;
subFolder_num       = 1;
baseFolder_name     = dir0(baseFolder_num).name;
dirBaseFolder       = dir(strcat(baseFolder_name,filesep,'1*_mat_Or'));
numSubFolders       = size(dirBaseFolder,1);
subFolder_name      = dirBaseFolder(subFolder_num).name;
caseName            = strcat(baseFolder_name,filesep,subFolder_name);             
disp(caseName)
caseNumber = 0;
% for counterCase = 1:size(DataARC_Datasets,1)
%         if strcmp(DataARC_Datasets{counterCase,2}(1:end-1),caseName(1:end-7))
%             caseNumber = counterCase;        
%         end
% end   

for counterCase = 1:size(DataARC_Datasets,1)
    % As the filesep can cause problems compare the base folder
    % first and then the filename. The filesep should be the
    % 11th e.g.  LSM_180907\180907_Timelapse_9
    
    if strncmp(DataARC_Datasets{counterCase,2}(1:end-1),caseName(1:end-7),10)
        if strcmp(DataARC_Datasets{counterCase,2}(12:end-1),caseName(12:end-7))
            caseNumber = counterCase;
        end
    end
end

if caseNumber~=0
  channelDistribution = [DataARC_Datasets{caseNumber,(12:19)}];
else
    clear channelDistribution;
end
%%
% Extract measurements 
% add the calibration parameters to obtain results in um/sec instead of
% just pixels/frame   [DataARC_Datasets{ca seNumber,22:23}]
counterCells =1;
[cell_metrics, nuclei_metrics] = measurementExtraction(caseName,counterCells,channelDistribution,1,[],[DataARC_Datasets{caseNumber,22:23}]);
%[cell_metrics, nuclei_metrics] = measurementExtraction(dir0,neutrophilNumber,ChannelDistribution,displayImages,selectedTimeFrame)
close all
%measurementExtraction('LSM_170403/170403_Timelapse_9_mat_Or',[],1);
%%
% quick partition
if exist('channelDistribution','var')
    quickPartition(caseName,channelDistribution,86)
else
    quickPartition(caseName)
end
%%



% quick ADD partition
if exist('channelDistribution','var')
    addPartition(caseName,channelDistribution,1)
else
    addPartition(caseName)
end



%% quick visualise

quickVisualise(caseName)



%%

figure
subplot(511)
plot([cell_metrics.numEndP],'k-','linewidth',2)
axis tight;grid on
title('Num Endpoints')

subplot(512)
plot([cell_metrics.NewArea]./[cell_metrics.Area],'b-','linewidth',2)
axis tight;grid on
title('New Area/Cell Area')

subplot(513)
plot([cell_metrics.MinorAxisLength]./[cell_metrics.MajorAxisLength],'m-','linewidth',2)
axis tight;grid on
title('Min Axis/Maj Axis')

subplot(514)
plot([nuclei_metrics.PositionS]/100,'r-','linewidth',2)
axis tight;grid on
title('Position Nuclei')


subplot(515)
plot([cell_metrics.Dist],'color',[0 0.6 0],'linewidth',2)
axis tight;grid on
title('Centroid Movement')
