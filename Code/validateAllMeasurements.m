%clear validateChannels DataARC_Datasets
clear all 
close all



%%
if strcmp(filesep,'/')
    % Running in Mac
    cd '/Users/ccr22/Academic/ARC_Grant/ARC'   % Run in PowerMac
    %cd /Volumes/TOSHIBA/ARC/                    % Run in Laptop
    load('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Datasets/DataARC_Datasets_2019_05_03.mat')
    %Output Folder
    dataOutFolder   = '/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Results/Metrics_2019_04_25';
    
else
    % running in windows
    cd ('D:\Acad\ARC\ARC')
    %load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2019_05_03.mat')
    load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2020_03_23.mat')
    %Output Folder
    %dataOutFolder   = 'D:\OneDrive - City, University of London\Acad\ARC_Grant\Results\Metrics_2019_04_25\';
    dataOutFolder   = 'D:\OneDrive - City, University of London\Acad\ARC_Grant\Results\Metrics_2020_03_23\';
end

%%
%load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\newCases_2020_03_23.mat')

if exist('DataARC_Datasets','var')
    % Only start if the data set info is available for the
    % ChannelDistribution and number of cells
    
    % Read all folders with data (grouped by date of acquisition)
    dir0            = dir ('LSM*');
    numFolders_0    = size(dir0,1);
    % There are files and folders, so select **manually** the important ones to validate
    % the data sets
    disp(' - - - - - - - - - - - - - - - - - - -')
    for baseFolder_num  =57:numFolders_0 % 5;
        % Read the timelapses in every day
        baseFolder_name = dir0(baseFolder_num).name;
        dirBaseFolder       = dir(strcat(baseFolder_name,filesep,'2*_mat_Or'));

        numSubFolders   = size(dirBaseFolder,1);
        
        for subFolder_num   =  1:numSubFolders  % 18;
            % Read every
            subFolder_name  = dirBaseFolder(subFolder_num).name;
            caseName            = strcat(baseFolder_name,filesep,subFolder_name);
            
            %disp(caseName)
            
            % Find case of all data sets (if DataARC_Datasets exists)
            
            caseNumber = 0;
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
                numCells = DataARC_Datasets{caseNumber,20};
            else
                clear channelDistribution;
            end
            %
            % Finally loop in case there are several cells in every set
            for counterCells =1:numCells
                disp(strcat(caseName,'  - ',num2str(counterCells)))
                clear cell_metrics nuclei_metrics
                % only process IF it is not discarded
                if (DataARC_Datasets{caseNumber,21}>0)
                    % add the calibration parameters to obtain results in um/sec instead of
                    % just pixels/frame   [DataARC_Datasets{caseNumber,22:23}]
                    [cell_metrics, nuclei_metrics] = measurementExtraction(caseName,counterCells,channelDistribution,1,[],[DataARC_Datasets{caseNumber,22:23}],dataOutFolder);
                    close all
                end
            end
            
            
            
            %numSet = numSet+1;
            
            
        end
    end
    
end




disp(' - - - - - - - - - - - - - - - - - - -')

%%
