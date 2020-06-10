function [cell_metrics, nuclei_metrics] = measurementExtraction(dir0,neutrophilNumber,ChannelDistribution,displayImages,selectedTimeFrame,calibrationParameters,dataOutFolder)

% The calibration parameters are used to obtain results in um/sec instead of
% just pixels/frame   [DataARC_Datasets{caseNumber,22:23}] and are
%    [um/pix  sec/frame]

%% Initial parsing of input data
if ~exist('dir0','var')
    disp('No input folder detected, discard case');
    return;
    %dir0       = 'LSM_170526/170526_timelapse17_mat_Or/';
end
%% read the folder where the data time frames are stored
dir1                        = dir(strcat(dir0,'/T*.mat'));
% number of time frames
numFrames                   = size(dir1,1);
% Read the first time frame to get the dimensions
load(strcat(dir0,'/',dir1(1).name));

% Calculate the distribution of the channels based on the dimensions
% Determine dimensions of the data sets
[rows,columns,numSlices ]   = size(dataIn);

if (~exist('ChannelDistribution','var'))||(isempty(ChannelDistribution))
    disp('---- No Channel Distribution detected ----')
    numSlicesPerCh              = numSlices/4;
    ChannelDistribution         = [1 numSlicesPerCh 1+numSlicesPerCh 2*numSlicesPerCh 1+2*numSlicesPerCh 3*numSlicesPerCh];
end

if ~exist('neutrophilNumber','var')
    neutrophilNumber       = 1;
end

%% Process the mask 
find_Or                         = strfind(dir0,'_Or');
dirMask                         = strcat(dir0(1:find_Or),'mask/');

if ~isdir(dirMask)
    % if there is no mask directory, create partitions with all the image
    % as a single partition
    %disp('No mask folder, discard case');
    %return;
    partitions                  = ones(rows,columns);
else
    % Check that the partition matches with the selected cell
    load(strcat(dirMask,'/',dir1(1).name)); 
end

highestCellNumber = max(partitions(:));
if neutrophilNumber > highestCellNumber
    partitionNotPresent         = 1;
    if partitionNotPresent ==1
        disp('Cell number does not contain a partition')
        cell_metrics            = [];
        nuclei_metrics          = [];
        return
    end 
end


if ~exist('displayImages','var')
    displayImages               = 1;
end

% Define time frame of experiment
if (~exist('selectedTimeFrame','var'))||(isempty(selectedTimeFrame))
    selectedTimeFrame = 1:numFrames;
end

gaussFilt = fspecial('Gaussian',[7 7],5);
%%
hFig=figure;
jet2=jet;jet2(1,:)=[0 0 0];colormap(jet2)
set(hFig,'position',[35         366        1391         432])
tic;


cell_metrics2(selectedTimeFrame(end)).NewArea               = NaN;
cell_metrics2(selectedTimeFrame(end)).OldArea               = NaN;
cell_metrics2(selectedTimeFrame(end)).Dist                  = NaN;
cell_metrics2(selectedTimeFrame(end)).Angle                 = NaN;
cell_metrics2(selectedTimeFrame(end)).Extent                = NaN;
cell_metrics2(selectedTimeFrame(end)).numBranches           = NaN;
cell_metrics2(selectedTimeFrame(end)).numEndP               = NaN;
cell_metrics2(selectedTimeFrame(end)).avTortuosity          = NaN;

nuclei_metrics2(selectedTimeFrame(end)).NewArea             = NaN;
nuclei_metrics2(selectedTimeFrame(end)).OldArea             = NaN;
nuclei_metrics2(selectedTimeFrame(end)).Dist                = NaN;
nuclei_metrics2(selectedTimeFrame(end)).Dist                = NaN;
nuclei_metrics2(selectedTimeFrame(end)).Angle               = NaN;
nuclei_metrics2(selectedTimeFrame(end)).Position            = NaN;
nuclei_metrics2(selectedTimeFrame(end)).PositionS           = NaN;
nuclei_metrics2(selectedTimeFrame(end)).numBranches_Thin    = NaN;
nuclei_metrics2(selectedTimeFrame(end)).numEndP_Thin        = NaN;
nuclei_metrics2(selectedTimeFrame(end)).numBranches_Skel    = NaN;
nuclei_metrics2(selectedTimeFrame(end)).numEndP_Skel        = NaN;
nuclei_metrics2(selectedTimeFrame(end)).avTortuosity        = NaN;


% Iterate over number of frames or selected time frame from input
for k = selectedTimeFrame   %  1:numFrames
    
    %%
    % k=1;
    % t2=toc;
    %disp([ k t2] )
    load(strcat(dir0,'/',dir1(k).name));
    if isdir(dirMask)
        % only read if there is a folder with masks
        load(strcat(dirMask,'/',dir1(k).name));
    end
    
    
    %%
    currentMask = (partitions==neutrophilNumber);
    if max(max(bwlabel((currentMask))))
        %currentMask = imclose(currentMask,strel('disk',60));
    end
    numPartitions       = max(partitions(:));
    
    if (neutrophilNumber > numPartitions)||(~ismember(neutrophilNumber,unique(partitions)) )
        %skip in case the neutrophil is not present in the current frame
        %
        cell_metrics2(k).NewArea        = NaN;
        cell_metrics2(k).OldArea        = NaN;
        cell_metrics2(k).Dist           = NaN;
        cell_metrics2(k).Angle          = NaN;
        cell_metrics2(k).Extent         = NaN;
        cell_metrics2(k).numBranches    = NaN;
        cell_metrics2(k).numEndP    	= NaN;
        cell_metrics2(k).avTortuosity   = NaN;
        cell_metrics2(k).touchBorder    = NaN;
        
        nuclei_metrics2(k).NewArea      = NaN;
        nuclei_metrics2(k).OldArea      = NaN;
        nuclei_metrics2(k).Dist         = NaN;
        nuclei_metrics2(k).Dist         = NaN;
        nuclei_metrics2(k).Angle        = NaN;
        nuclei_metrics2(k).Position     = NaN;
        nuclei_metrics2(k).PositionS    = NaN;
        
        
    else
        
        %         The channel assignation has been rather inconsistent:
        %             (Using the example of 170412_timelapse8 to check)
        %             files with 4 channels are:
        %             C0 are green
        %             C1 blue
        %             C3 red
        %             C4 DIC
        %
        %             For THREE channel data
        %             C0 green
        %             C1 red
        %             C2 DIC
        %
        %             Except for THREE channel data from 170727 (an exception to the rule, remembering in this case that blue and green colours were not properly separated in acquisition anyway)
        %             C0 blue
        %             C1 green
        %             C2 red
        %             No DIC
        %
        %       Thus the matrix DataARC_Datasets will have the following arrangment:
        %               (:,12:13) g_i   g_f
        %               (:,14:15) b_i   b_f
        %               (:,16:17) r_i   r_f
        %               (:,18:19) DIC_i DIC_f
        %       If no channels are present, the field will be empty
        %
        
        
        
        % get the mean value of each channel
        if (ChannelDistribution(1)~=0)
            g_Channel           = dataIn(:,:,ChannelDistribution(1):ChannelDistribution(2));
            g0                  = imfilter(mean(g_Channel,3),gaussFilt,'replicate');
            g                   = (currentMask).*g0;
        else
            g_Channel       = zeros(rows,columns,1+ChannelDistribution(1)-ChannelDistribution(2));
            g               = zeros(rows,columns);
            g0              = g;
        end
        %b = mean(dataIn(:,:,ChannelDistribution(3):ChannelDistribution(5)),3);
        %r = mean(dataIn(:,:,ChannelDistribution(5):ChannelDistribution(6)),3);
        
        %if numel(ChannelDistribution)>2
        if (ChannelDistribution(3)~=0)
            b_Channel       = dataIn(:,:,ChannelDistribution(3):ChannelDistribution(4));
            b0              = imfilter(mean(b_Channel,3),gaussFilt,'replicate');
            b               = ((currentMask).*b0);
        else
            b_Channel       = zeros(rows,columns,1+ChannelDistribution(4)-ChannelDistribution(3));
            b               = zeros(rows,columns);
            b0              = b;
        end
        %if numel(ChannelDistribution)>4
        if (ChannelDistribution(5)~=0)
            r_Channel       = dataIn(:,:,ChannelDistribution(5):ChannelDistribution(6));
            r0              = imfilter(mean(r_Channel,3),gaussFilt,'replicate');
            r               = ((currentMask).*r0);
        else
            r_Channel       = zeros(rows,columns,1+ChannelDistribution(6)-ChannelDistribution(5));
            r               = zeros(rows,columns);
            r0              = r;
        end
        %imagesc(composite(:,:,:,k))
        
        % the maxima are only those in the regions of the cells, this avoids those rare
        % blue regions in some images
        max_R               = max(max(r0.*(partitions>0)));
        max_G               = max(max(g0.*(partitions>0)));
        max_B               = max(max(b0.*(partitions>0)));
        
        
        r1                  = r/max_R;
        g1                  = g/max_G;
        b1                  = b/max_B;
        composite(:,:,1)    = r1 ;
        composite(:,:,2)    = g1;
        composite(:,:,3)    = b1;
        composite0(:,:,1)   = r0/max_R;
        composite0(:,:,2)   = g0/max_G;
        composite0(:,:,3)   = b0/max_B;
        
        thres_R             = graythresh(r(currentMask)/max_R);
        thres_G             = graythresh(g(currentMask)/max_G);
        thres_B             = graythresh(b(currentMask)/max_B);
        % For normal cases use 0.7
       % r2                  = r1 > max(0.10,(0.7*thres_R));
       
        % For the tethered cases use lower levels
        r2                  = r1 > max(0.06,(0.1*thres_R));
        g2                  = g1 > max(0.2,(0.9*thres_G));
        b2                  = b1 > max(0.2,(0.9*thres_B));
       
        
        %thresCell           = max(0.2,min(0.9*[thres_R thres_G thres_B]));
        maxIntProjCell      = max(composite,[],3);
        maxIntProjNucl      = max(composite(:,:,2:3),[],3);
        
        % Find the Whole Cell by combining the three channels and thresholding together
        %cellMask            = imclose(maxIntProjCell>thresCell,strel('disk',6));
        
        cellMask_init            = imclose((r2+b2+g2)>0,strel('disk',6));
        
        % there should only be one cell so discard small regions, keep the largest one
        cellMask_L          = bwlabel(cellMask_init);
        cellMask_R          = regionprops(cellMask_L,maxIntProjCell,'Area','MeanIntensity');
        [areas_1,areas_2]   = sort([cellMask_R.Area]);
        [inten_1,inten_2]   = sort([cellMask_R.MeanIntensity]);
        [comb_1,comb_2]     = sort([cellMask_R.Area].*[cellMask_R.MeanIntensity]);
        try
            % As the cells have some thin pseudopods, the masks of partitions
            % should be taken into account to prevent the joining of regions
            cellMask            = currentMask.*imdilate(ismember(cellMask_L,comb_2(end)),ones(3));
            %cellMask            = (partitions>0).*imdilate(ismember(cellMask_L,comb_2(end)),ones(3));
            
            %ensure there are no holes or spurious in the cell Mask
            cellMask            = imfill(cellMask,'holes');
            cellMask            = bwmorph(cellMask,'spur',1);
            
            
        catch
            qqq=1;
        end
        [cellMask_L2,num_CM]          = bwlabel(cellMask);
        if num_CM>1
            cellMask_R2          = regionprops(cellMask_L2,'Area');
            cellMask            = ismember(cellMask_L2,find([cellMask_R2.Area]== max([cellMask_R2.Area])));
        end
        cellMaskBor         = zerocross(cellMask);
        
        %cellMask_F2          = ismember(cellMask_L,areas_2(end));
        %thresNucl           = max(0.3,1.1*graythresh(maxIntProjNucl(cellMask_F)));
        
        % Repeat the process for the nuclei with G+B remove areas where red is very
        % strong
        g3                  = g1 > max(0.2,1.2*(thres_G));
        b3                  = b1 > max(0.2,1.2*(thres_B));
        %b3                  = b1 > max(0.2,0.8*(thres_B));
        r4                  = r1 > min(0.5,(3*thres_R));
        
        r5                  = r1.*cellMask;
        g5                  = g1.*cellMask;
        b5                  = b1.*cellMask;
        
        
        
        % the Nuclei is defined HERE, in some cases, it can be split so
        % check the conditions. 2019_04_08  Nuclei is sometimes outside the
        % cell, change so that it is ALWAYS inside
        
        % This only takes into account the blue/green channels
        %nucleiMask          = imopen( (b3|g3) ,strel('disk',2));
        nucleiMask          = cellMask.* imopen( (b3|g3) ,strel('disk',2));
        
        
        % ensure that the mask does not have holes
        nucleiMask          = imfill( nucleiMask,'holes');
        
        % This takes into account the red channel and suppresses the
        % intensity in that area.
        %nucleiMask          = imopen( (b3|g3).*((max(b5,g5))>r5) ,strel('disk',2));
        
        % IF there is more than one region of nuclei ... keep the largest?
        [nucleiMaskL,numNuclParts]= bwlabel (nucleiMask);
        
        if numNuclParts>1
            % If there is more than one region of Nuclei analyse and decide what to
            % do:
            % This rule select the LARGEST section of the nuclei, works well when
            % there is one region considerably larger (e.g. 540 . 120) than another,
            % however, if the difference is smaller (312 v 310) it may select the
            % incorrect one. Indeed, these should not be separated at all.
            nucleiMaskP     = regionprops(nucleiMaskL,'area');
            [~,sizeNucPart] = sort([nucleiMaskP.Area],'descend');
            nucleiMask      = (nucleiMaskL==sizeNucPart(1));
            
            % Consider not only size, but the intensity of the region
            if all(isnan(b1(:)))
                combIntensity = g1;
            elseif all(isnan(g1(:)))
                combIntensity = b1;
            else
                combIntensity = g1+b1;
            end
            nucleiMaskP     = regionprops(nucleiMaskL,combIntensity,'area','meanintensity','maxintensity','centroid');
            % If size is similar, or intensity is similar, then these should not be
            % split but a single region,
            % Find the relative sizes and intensities and combine in a single
            % variable of importance. relImportance ==1 means it is the largest and
            % brightest
            relArea         = (([nucleiMaskP.Area]/max([nucleiMaskP.Area])))';
            relIntens       = (([nucleiMaskP.MeanIntensity]/max([nucleiMaskP.MeanIntensity])))';
            relImportance   = (relArea+relIntens)/2;
            
            % Remove sections with relImportance less than 0.5, anything that is not
            % at least the brightest or the largest
            thresRelImportance = 0.6;
            if (sum(relImportance>thresRelImportance)==1)
                % only one important region, keep
                % beep
                nucleiMask      = (nucleiMaskL==find(relImportance>thresRelImportance));
            else
                % There are at least 2 regions with relative importance > 0.5
                % merge them with a close, hopefully these will get connected
                beep
                nucleiMask      = imclose( ismember(nucleiMaskL,find(relImportance>thresRelImportance))  ,strel('disk',20));
                
                if sum(sum((nucleiMaskL>0)~=nucleiMask))
                    % the closing does not connect sections, these must be appart
                    % thus discard and keep the one with highest relImportance
                    [~,importNucPart] = sort(relImportance,'descend');
                    nucleiMask      = (nucleiMaskL==importNucPart(1));
                    
                    %nucleiMask      = (nucleiMaskL==find(relImportance>thresRelImportance));
                end
                
                disp(strcat('-------',dir0))
                %                             figure(2);imagesc(nucleiMaskL+nucleiMask)
                %                             figure(3);imagesc(composite)
                %                             figure(1)
                %                                                         q=1;
            end
            
            
            
            
            
            clear nucleiMaskP sizeNucPart nucleiMaskL numNuclParts;
        end
        
        % Obtain the perimeter of the Nuclei
        nucleiMaskBor       = zerocross(nucleiMask);
        
        % nucleiMask          = imopen((b3|g3).*(1-r4),strel('disk',2));
        
        
        % there should only be one cell so discard small regions, keep the largest one
        [nucleiMask_L,numRegNuc]          = bwlabel(nucleiMask);
        if numRegNuc>1
            nucleiMask_R          = regionprops(nucleiMask_L,max(composite(:,:,2:3),[],3),'Area','MeanIntensity');
            [areas_1,areas_2]   = sort([nucleiMask_R.Area]);
            [inten_1,inten_2]   = sort([nucleiMask_R.MeanIntensity]);
            [comb_1,comb_2]       = sort([nucleiMask_R.Area].*[nucleiMask_R.MeanIntensity]);
            nucleiMask_F          = ismember(nucleiMask_L,comb_2(end));
            nucleiMask = nucleiMask_F;
        end
        
        
        %  imagesc(composite)
        %[thres_R thres_G thres_B ]
        %save(strcat(dirMask,dir1(k).name),'partitions');
        
        
        cellMask_t(:,:,k)       = cellMask;
        nucleiMask_t(:,:,k)     = nucleiMask;
        
        % Measurements to be extracted:
        % Cell
        % 1 Area/volume 2 elongation 3 max/min axis 4 speed 5 direction, 6 branches
        
        cell_metrics(k)     = regionprops(cellMask,'Area','Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation','BoundingBox','Extrema');
        cellSkel            = bwmorph(cellMask,'thin','inf');
        
        % Nuclei
        % 1 Area/volume 2 elongation 3 max/min axis 4 speed 5 direction, 6 branches
        
        % If nuclei exist, define skeletons that do not overlap with nuclei
        % and its boundary
        if max(nucleiMask(:))>0
            try
            nuclei_metrics(k)   = regionprops(nucleiMask,'Area','Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation');
            catch
                qqq=1;
            end
            nucleiSkel          = bwmorph(nucleiMask,'thin','inf');
            nucleiSkel3         = bwmorph(nucleiMask,'skel','inf');
            cellSkel2           = bwlabel(cellSkel.*(1-nucleiMask));
            cellSkel3           = cellSkel.*(1-zerocross(nucleiMask)).*(1-nucleiSkel);
            nucleiSkel2         = nucleiSkel.*(1-cellSkel);
            nucSkelEndPoints       = bwmorph(nucleiSkel,'endpoints');
            nucSkelEndPoints3      = bwmorph(nucleiSkel3,'endpoints');
            nucSkelBranchPoints    = bwmorph(nucleiSkel,'branchpoints');
            nucSkelBranchPoints3   = bwmorph(nucleiSkel3,'branchpoints');
            
            numEndpointsNucSkel    = sum(sum(nucSkelEndPoints));
            numEndpointsNucSkel3   = sum(sum(nucSkelEndPoints3));
            numBranchpointsNucSkel    = sum(sum(nucSkelBranchPoints));
            numBranchpointsNucSkel3   = sum(sum(nucSkelBranchPoints3));
            
            
            numPixNucSkel          = sum(sum(nucleiSkel));
            
            % Calculate a broken skeleton with only the segments between branch points
            nucleiSkelPerSegment      = nucleiSkel .*(1-imdilate(nucSkelBranchPoints,ones(3)));
            nucleiSkelPerSegment_L    = bwlabel( nucleiSkelPerSegment);
            nucleiSkelPerSegment_P    = regionprops( nucleiSkelPerSegment_L,'Area');
            % Calculate the "forkness" of the skeleton, ratio of the longest segment over
            % the total. A single branch will be 1, whilst a skeleton with three
            % segments, more or less equal will be 1/3
            nucleiForkness            = max([nucleiSkelPerSegment_P.Area])/sum([nucleiSkelPerSegment_P.Area]);
            
            
            
            
        else
            cellSkel2           = cellSkel;
            cellSkel3           = cellSkel;
            nucleiSkel          = zeros(rows, columns);
            nucleiSkel3         = zeros(rows, columns);
            numEndpointsNucSkel =nan;
            numBranchpointsNucSkel   =nan;
            numEndpointsNucSkel3 =nan;
             numBranchpointsNucSkel3   =nan;
             nucleiForkness            = nan;

            
        end
        numSegmentsCellSkel     = max(cellSkel2(:));
        cellSkelEndPoints       = bwmorph(cellSkel,'endpoints');
        cellSkelBranchPoints    = bwmorph(cellSkel,'branchpoints');
        numEndpointsCellSkel    = sum(sum(cellSkelEndPoints));
        numPixCellSkel          = sum(sum(cellSkel));
        % Find the locations of all ENDPOINTS of the skeleton
        [rEndPoints,cEndPoints] = find(cellSkelEndPoints);
        try
            rOrigin                 = rEndPoints(1);
        catch
            qqq=1;
        end
        cOrigin                 = cEndPoints(1);
        %[rOrigin,cOrigin]       = find(cellSkelEndPoints,1);
        
        % Calculate the average distance between endpoints
        % 1->2, 2->3, ...n-1->n,n->1
        endPointPositions       = [[rEndPoints cEndPoints];[rEndPoints(1) cEndPoints(1)]];
        AvDistEndPoints         = sum(sqrt(sum((diff(endPointPositions)).^2,2)));
        propsCellSkel           = regionprops(cellSkel,'perimeter');
        AvPerimeterCellSkel     = propsCellSkel.Perimeter;
        propsNucSkel            = regionprops(nucleiSkel,'perimeter');
        if max(nucleiSkel(:))==0
            AvPerimeterNucSkel      = 0;
        else
            AvPerimeterNucSkel      = propsNucSkel.Perimeter;
        end
        % Calculate a broken skeleton with only the segments between branch points
        cellSkelPerSegment      = cellSkel .*(1-imdilate(cellSkelBranchPoints,ones(3)));
        cellSkelPerSegment_L    = bwlabel( cellSkelPerSegment);
        cellSkelPerSegment_P    = regionprops( cellSkelPerSegment_L,'Area');
        % Calculate the "forkness" of the skeleton, ratio of the longest segment over
        % the total. A single branch will be 1, whilst a skeleton with three
        % segments, more or less equal will be 1/3
        cellForkness            = max([cellSkelPerSegment_P.Area])/sum([cellSkelPerSegment_P.Area]);
        
        
        
        % Calculate distance along the skeleton of the Cell, not just BWDIST as
        % there may be branches that should have the same distance along
        % themselves, also a curved line should have long distances even if it is
        % close Euclideanly from the origin.
        
        % Use bwtraceboundary
        clockWise_Boun     = bwtraceboundary(cellSkel,[rOrigin,cOrigin],'e',8,Inf,'clockwise');
        counterCW_Boun     = bwtraceboundary(cellSkel,[rOrigin,cOrigin],'e',8,Inf,'counterclockwise');
        
        dist_CW = (1+numPixCellSkel*2)*zeros(size(cellSkel));
        dist_CCW = dist_CW;
        
        % traverse the two
        for counter_Boun = 1:size(clockWise_Boun,1)
            dist_CW (clockWise_Boun(counter_Boun,1),clockWise_Boun(counter_Boun,2)) = numPixCellSkel*2-counter_Boun;
            dist_CCW(counterCW_Boun(counter_Boun,1),counterCW_Boun(counter_Boun,2)) = numPixCellSkel*2-counter_Boun;
        end
        dist_CellSkel = min(dist_CW,dist_CCW);
        
        % Calculate alignment of skeletons, distance between them, their lengths and
        % ratios
        
        cellSkel_dist                   = bwdist(cellSkel);
        skelAlignment                   = double(sum(sum(cellSkel_dist.*nucleiSkel)));
        
        
        
        
        cell_metrics2(k).numBranches    = numSegmentsCellSkel;
        cell_metrics2(k).numEndP        = numEndpointsCellSkel;
        cell_metrics2(k).skelPerim      = AvPerimeterCellSkel;
        cell_metrics2(k).AvDistEndP     = AvDistEndPoints;
        cell_metrics2(k).avTortuosity   = AvDistEndPoints/AvPerimeterCellSkel;
        cell_metrics2(k).forkness       = cellForkness;
        cell_metrics2(k).skelAlignment  = skelAlignment;
        %detect if the cell is touching the borders
        if (any(cell_metrics(k).Extrema(:)<=1)) || (any(cell_metrics(k).Extrema(:,2)>=rows))  || (any(cell_metrics(k).Extrema(:,1)>=columns))
            cell_metrics2(k).touchBorder    = 1;
        else
            cell_metrics2(k).touchBorder    = 0;
        end
        
        % Calculate the average position of the nuclei inside the cell
        
        length_CellSkel     = max(dist_CellSkel(:));
        
        positionNuclei      =(nucleiMask.* dist_CellSkel);
        positionNuclei(positionNuclei==0)=[];
        % scale to 0-100 to know where it is
        %positionNucleiR = round(100*mean(positionNuclei)/length_CellSkel);
        % scale to [-1,+1] to know where it is
        positionNucleiR = (-1+2*mean(positionNuclei)/length_CellSkel);
        
        nuclei_metrics2(k).PositionS = positionNucleiR;
        
        nuclei_metrics2(k).forkness     = nucleiForkness;
        nuclei_metrics2(k).skelPerim    = AvPerimeterNucSkel;
        
        
        
        % Combined
        %    position of the nuclei in the cell
        %   velocity of the cell as centroid(t+1) - centroid (t)
        %   movement of the cell as mask(t+1) ~= mask(t)
        if k>1
            if ~isempty(cell_metrics(k-1).Centroid)
                try
                    movementVect         = cell_metrics(k).Centroid -cell_metrics(k-1).Centroid;
                catch
                    qqq=1;
                end
                movementVect2            = movementVect(1) +1i*movementVect(2);
                movementAreaC            = 2*cellMask_t(:,:,k)+cellMask_t(:,:,k-1);
                distFromPrevious         = bwdist(imopen(movementAreaC==1,ones(2))).*cellMask_t(:,:,k);
                cell_metrics2(k).NewArea = sum(sum(  movementAreaC==2 ));
                cell_metrics2(k).OldArea = sum(sum(  movementAreaC==1 ));
                cell_metrics2(k).Dist    = abs(movementVect2);   %sqrt(sum((movementVect).^2));
                cell_metrics2(k).Angle   = angle(movementVect2);
                cell_metrics2(k).Extent  = max(distFromPrevious(:));
                
                
                if max(nucleiMask(:)>0)
                    if ~isempty(nuclei_metrics(k-1).Centroid)
                        try
                            movementVect             = nuclei_metrics(k).Centroid -nuclei_metrics(k-1).Centroid;
                        catch
                            qqq=1;
                        end
                        movementVect2            = movementVect(1) +1i*movementVect(2);
                        movementAreaN             = 2*nucleiMask_t(:,:,k)+nucleiMask_t(:,:,k-1);
                        nuclei_metrics2(k).NewArea = sum(sum(  movementAreaN==2 ));
                        nuclei_metrics2(k).OldArea = sum(sum(  movementAreaN==1 ));
                        nuclei_metrics2(k).Dist    = sqrt(sum((movementVect).^2));
                        nuclei_metrics2(k).Dist    = abs(movementVect2);   %sqrt(sum((movementVect).^2));
                        nuclei_metrics2(k).Angle   = angle(movementVect2);
                        % this is trying to calculate the relative position of the
                        % nuclei with respect to the cell by calculating the
                        % difference between the current cell and the previous one,
                        % this works well ONLY if the cell is advancing forward
                        % relatively fast so that the difference between positions is
                        % clearly in front and back. Does not work when there are
                        % small movements at the sides.
                        nucleiPosition0      =  distFromPrevious.*(nucleiMask_t(:,:,k));
                        nucleiPosition1      = nucleiPosition0(:);
                        nucleiPosition1(nucleiPosition1==0)=[];
                        nucleiPosition       = mean(nucleiPosition1)/cell_metrics2(k).Extent;
                        nuclei_metrics2(k).Position = nucleiPosition;
                        
                        
                        nuclei_metrics2(k).numBranches_Thin = numBranchpointsNucSkel;
                        nuclei_metrics2(k).numEndP_Thin     = numEndpointsNucSkel;
                        nuclei_metrics2(k).numBranches_Skel = numBranchpointsNucSkel3;
                        nuclei_metrics2(k).numEndP_Skel     = numEndpointsNucSkel3;

                        
                        
                        
                    end
                end
                
                rowAx=[-5+cell_metrics(k).BoundingBox(2) 5+cell_metrics(k).BoundingBox(2)+cell_metrics(k).BoundingBox(4)];
                colAx=[-5+cell_metrics(k).BoundingBox(1) 5+cell_metrics(k).BoundingBox(1)+cell_metrics(k).BoundingBox(3)];
                
                
                
                %figure(2)
                %    imagesc(composite);%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
                
                
                
                
                maxFontSize =9;
                if displayImages==1
                    % Display section
                    subplot(141)
                    imagesc(imdilate(zerocross(cellMask-nucleiMask),ones(1))+composite0);%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
                    %figure(3)
                    % imagesc(imdilate(zerocross(partitions),ones(3))+composite0);%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
                    %figure(1)
                    %  imagesc(();%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
                    %imagesc(imdilate(zerocross(cellMask),ones(3))+max(composite,[],3));colorbar
                    title(strcat(dir0,'   frame=',num2str(k),'/',num2str(numFrames),';  cell=',num2str(neutrophilNumber),'/',num2str(numPartitions)),'interpreter','none','fontsize',maxFontSize)
                    drawnow
                    pause(0.005)
                    
                    subplot(142)
                    imagesc(movementAreaC)
                    axis([ colAx rowAx])
                    
                    title('Cyan (t-1), Yellow(t), Brown (t,t-1)','interpreter','none','fontsize',maxFontSize)
                    
                    %subplot(142)
                    % Display the distance from the previous region as compared with the
                    % current nuclei and cell mask
                    %             imagesc(distFromPrevious.*(cellMask_t(:,:,k)-nucleiMask_t(:,:,k)))
                    %             axis([ colAx rowAx])
                    %             colorbar
                    %             title('Distance from area in (t-1)')
                    %
                    subplot(143)
                    %imagesc(cellMask+cellSkel+2*nucleiMask)
                    imagesc(4*nucleiSkel  +imdilate(cellSkelEndPoints,ones(3)) +cellSkel3 + (1+numSegmentsCellSkel)*(cellMaskBor)+(2+numSegmentsCellSkel)*(nucleiMaskBor))
                    axis([ colAx rowAx])
                    title(strcat('Num Branches Cell Skeleton = ',num2str(numEndpointsCellSkel)),'fontsize',maxFontSize)
                    
                    
                    %                     subplot(144)
                    %                     imagesc((1-nucleiMaskBor).*dist_CellSkel+ (length_CellSkel)*(cellMaskBor)+(numPixCellSkel/2)*(nucleiMaskBor))
                    %                     axis([ colAx rowAx])
                    %                     title(strcat('Rel Position Nuc = ',num2str(positionNucleiR),'/100'))
                    
                    
                    r6                  = zeros(rows,columns,3);
                    r6(:,:,1)           = r5+cellMaskBor+nucleiMaskBor+nucleiSkel;
                    r6(:,:,2)            = cellMaskBor+nucleiMaskBor+nucleiSkel;
                    r6(:,:,3)            = cellMaskBor+nucleiMaskBor+nucleiSkel;
                    g6                  = zeros(rows,columns,3);
                    g6(:,:,2)           = g5+cellMaskBor+nucleiMaskBor+nucleiSkel3;
                    g6(:,:,1)            = cellMaskBor+nucleiMaskBor+nucleiSkel3;
                    g6(:,:,3)            = b5+ cellMaskBor+nucleiMaskBor+nucleiSkel3;
                    
                    r6(r6>1) =2;
                    g6(g6>1) =1;
                    subplot(244)
                    imagesc(r6)
                    title(strcat('End/Branch points (thin)= ',num2str(numEndpointsNucSkel),'/',num2str(numBranchpointsNucSkel)))
                    
                    axis([ colAx rowAx])
                    %
                    
                    subplot(248)
                    imagesc(g6)
                    title(strcat('End/Branch points (skel)= ',num2str(numEndpointsNucSkel3),'/',num2str(numBranchpointsNucSkel3)))
                    axis([ colAx rowAx])
                    F(k) = getframe(hFig);
                end
            else
                subplot(141)
                imagesc(imdilate(zerocross(cellMask-nucleiMask),ones(1))+composite0);%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
                
            end
        else
            subplot(141)
            imagesc(imdilate(zerocross(cellMask-nucleiMask),ones(1))+composite0);%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
            %figure(3)
            % imagesc(imdilate(zerocross(partitions),ones(3))+composite0);%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
            %figure(1)
            %  imagesc(();%+repmat(imdilate(cellBound==0,ones(3)),[1 1 3]))
            %imagesc(imdilate(zerocross(cellMask),ones(3))+max(composite,[],3));colorbar
            title(strcat(dir0,'   frame=',num2str(k),'/',num2str(numFrames),';  cell=',num2str(neutrophilNumber),'/',num2str(numPartitions)),'interpreter','none')
            drawnow
            F(k) = getframe(hFig);
            %colormap jet
            %             figure(2)
            %             imagesc(cellSkel)
            %             axis([ colAx rowAx])
        end
    end
    
    if k==85
        qqq=1;
    end
    
    qqq=1;
end

qqq=2;
%%
cell_metrics2(1).NewArea =0;
for k = selectedTimeFrame   %  1:numFrames
    cell_metrics(k).NewArea     = cell_metrics2(k).NewArea;
    cell_metrics(k).OldArea     = cell_metrics2(k).OldArea;
    cell_metrics(k).Dist        = cell_metrics2(k).Dist;
    cell_metrics(k).Angle       = cell_metrics2(k).Angle;
    cell_metrics(k).Extent      = cell_metrics2(k).Extent;
    cell_metrics(k).numBranches = cell_metrics2(k).numBranches;
    cell_metrics(k).numEndP     = cell_metrics2(k).numEndP;
    cell_metrics(k).avTortuosity= cell_metrics2(k).avTortuosity;
    cell_metrics(k).touchBorder = cell_metrics2(k).touchBorder;
    cell_metrics(k).rows        = rows;
    cell_metrics(k).columns     = columns;
    cell_metrics(k).numSlices   = numSlices;
    cell_metrics(k).forkness    = cell_metrics2(k).forkness;
    cell_metrics(k).skelAlignment = cell_metrics2(k).skelAlignment;
    cell_metrics(k).skelPerim   = cell_metrics2(k).skelPerim;
    cell_metrics(k).ratioSkels  = nuclei_metrics2(k).skelPerim/cell_metrics2(k).skelPerim;
    
    nuclei_metrics(k).NewArea   = nuclei_metrics2(k).NewArea;
    nuclei_metrics(k).OldArea   = nuclei_metrics2(k).OldArea;
    nuclei_metrics(k).Dist      = nuclei_metrics2(k).Dist;
    nuclei_metrics(k).Angle     = nuclei_metrics2(k).Angle;
    nuclei_metrics(k).Position  = nuclei_metrics2(k).Position;
    nuclei_metrics(k).PositionS  = nuclei_metrics2(k).PositionS;
    
    
    nuclei_metrics(k).numBranches_Thin  = nuclei_metrics2(k).numBranches_Thin;
    nuclei_metrics(k).numEndP_Thin      = nuclei_metrics2(k).numEndP_Thin    ;
    nuclei_metrics(k).numBranches_Skel  = nuclei_metrics2(k).numBranches_Skel;
    nuclei_metrics(k).numEndP_Skel      = nuclei_metrics2(k).numEndP_Skel    ;
    nuclei_metrics(k).forkness          = nuclei_metrics2(k).forkness    ;
    nuclei_metrics(k).skelPerim         = nuclei_metrics2(k).skelPerim;
    

    
end

qqq=2;


fields_cell = fieldnames(cell_metrics);
fields_nuc  = fieldnames(nuclei_metrics);
% change all empty fields to NaN
for k = selectedTimeFrame
    for counterF = 1:numel(fields_nuc)
        if isempty(getfield(nuclei_metrics(k),fields_nuc{counterF}))
            nuclei_metrics(k) = setfield(nuclei_metrics(k),fields_nuc{counterF},NaN);
        end
    end
    for counterF = 1:numel(fields_cell)
        if isempty(getfield(cell_metrics(k),fields_cell{counterF}))
            cell_metrics(k) = setfield(cell_metrics(k),fields_cell{counterF},NaN);
        end
    end
    
end


%% Calibrate to um/sec if there are calibrationParamters

if exist('calibrationParameters','var')
    % calibrationParameters is a 1 x 2 matrix with the values of
    %  [um/pix sec/frame]
    
    for k = selectedTimeFrame   %  1:numFrames
        %                                       pix                             um                              frames
        %                                       -----                           -----                           -----
        %                                       frame                           pix                             sec
        cell_metrics(k).Dist_um_s    = cell_metrics(k).Dist             *    calibrationParameters(1)/calibrationParameters(2) ;
        nuclei_metrics(k).Dist_um_s  = nuclei_metrics(k).Dist           *    calibrationParameters(1)/calibrationParameters(2) ;
        
        cell_metrics(k).Area_um_2    = cell_metrics(k).Area             *    calibrationParameters(1)*calibrationParameters(1) ;
        cell_metrics(k).MajAxis_um   = cell_metrics(k).MajorAxisLength  *    calibrationParameters(1) ;
        cell_metrics(k).MinAxis_um   = cell_metrics(k).MinorAxisLength  *    calibrationParameters(1) ;
        cell_metrics(k).Min_MajAxis  = cell_metrics(k).MinorAxisLength/cell_metrics(k).MajorAxisLength;
        
        nuclei_metrics(k).Area_um_2  = nuclei_metrics(k).Area           *    calibrationParameters(1)*calibrationParameters(1) ;
        nuclei_metrics(k).MajAxis_um = nuclei_metrics(k).MajorAxisLength*    calibrationParameters(1) ;
        nuclei_metrics(k).MinAxis_um = nuclei_metrics(k).MinorAxisLength*    calibrationParameters(1) ;
        nuclei_metrics(k).Min_MajAxis= nuclei_metrics(k).MinorAxisLength/nuclei_metrics(k).MajorAxisLength;
        
    end
    
end
%% Recalculate the position of the nuclei once all the path has been calculated

% This does not work when there is a NaN as centroid as it only passes one argument
% instead of 2 (x,y) so do a loop instead
% % Take all movement to calculate the route, use centroids
% cent1                   = [cell_metrics.Centroid];
% centx                   = cent1(1:2:end);
% centy                   = cent1(2:2:end);
for counterT = selectedTimeFrame
    if isnan(cell_metrics(counterT).Centroid)
        centx(counterT)     = NaN;
        centy(counterT)     = NaN;
    else
        try
        centx(counterT)     = cell_metrics(counterT).Centroid(1);
        centy(counterT)     = cell_metrics(counterT).Centroid(2);
        catch
            q=1;
        end
    end
end
% filter the positions to remove small variations that have a strong variation in
% direction
centy2                  = imfilter(centy,[1 1 1 1 1]/5,'replicate');
centx2                  = imfilter(centx,[1 1 1 1 1]/5,'replicate');
% calculate the movement and the angle of direction
movementVectC           = centx2(2:end) -centx2(1:end-1) +1i*(centy2(2:end) -centy2(1:end-1));
angleR                  = 180*angle(movementVectC)/pi;



% Using the angle, rotate the cell so that the direction of movement is always
% towards the right.
relPositionNuclei(numFrames)=0;
nuclei_metrics(1).PositionR = 0;
nuclei_metrics(1).angleR    = 0;

for k=2:numFrames
    if isnan(angleR(k-1))
        nuclei_metrics(k).PositionR  = NaN;
        nuclei_metrics(k).angleR    = NaN;
    else
        cellMaskR               = imrotate(cellMask_t(:,:,k),angleR(k-1));
        cellMaskR2              = max(cellMaskR,[],1);
        nucleiMaskR             = imrotate(nucleiMask_t(:,:,k),angleR(k-1));
        nucleiMaskR2            = max(nucleiMaskR,[],1);
        nucleiPositionR         = (nucleiMaskR2.*(cumsum(cellMaskR2)/sum(cellMaskR2)));
        nucleiPositionR(nucleiPositionR==0)=[];
        relPositionNuclei(k)       = 2*mean(nucleiPositionR)-1;
        %     subplot(2,10,k-60)
        %     imagesc((cellMask_t(:,:,k)+nucleiMask_t(:,:,k)))
        %     axis off
        %     subplot(2,10,k+10-60)
        %     imagesc(cellMaskR+nucleiMaskR)
        %     axis off
        %     drawnow;
        nuclei_metrics(k).PositionR  = relPositionNuclei(k);
        nuclei_metrics(k).angleR    = angleR(k-1);
    end
end
%% Calculate tortuosity of the PATH of the cell based on the positions

pathy=centy2([1:end]);
pathx=centx2([1:end]);
pathy(isnan(pathy))=[];
pathx(isnan(pathx))=[];

distanceBetweenExtremes         = sqrt((diff(pathy([1 end]))).^2+(diff(pathx([1 end]))).^2);

distanceOfPath                  = sum(sqrt((diff(pathy)).^2+(diff(pathx)).^2));


pathTortuosity                  = distanceOfPath/distanceBetweenExtremes;
cell_metrics(1).pathTortuosity  = pathTortuosity;
%% Save the results of cell/nuclei metrics


if ~exist('dataOutFolder','var')    
    MetricsOutName      = strcat(dir0(1:end-6),'cell_',num2str(neutrophilNumber),'_metrics.mat');
    CellMetricsOutName  = strcat(dir0(1:end-6),'cell_',num2str(neutrophilNumber),'_metrics.xls');
    figOutName = strcat(dir0(1:end-6),'cell_',num2str(neutrophilNumber),'_stats.jpg');
    videoOutName = strcat(dir0(1:end-6),'cell_',num2str(neutrophilNumber),'.mp4');
else
    MetricsOutName      = strcat(dataOutFolder,'metrics',filesep,dir0(12:end-6),'cell_',num2str(neutrophilNumber),'_metrics.mat');
    CellMetricsOutName  = strcat(dataOutFolder,'metrics_xls',filesep,dir0(12:end-6),'cell_',num2str(neutrophilNumber),'_metrics.xls');
    figOutName = strcat(dataOutFolder,'stats/',filesep,dir0(11:end-6),'cell_',num2str(neutrophilNumber),'_stats.jpg');
    videoOutName = strcat(dataOutFolder,'videos/',filesep,dir0(11:end-6),'cell_',num2str(neutrophilNumber),'.mp4');

%    dataOutFolder = '';
end


save(MetricsOutName,'cell_metrics','nuclei_metrics')
%writetable(struct2table(cell_metrics),CellMetricsOutName,'Sheet','Cell');
%writetable(struct2table(nuclei_metrics),CellMetricsOutName,'Sheet','Nuclei');


%%

if (displayImages==1)
    %% Print a picture and Generate a video
    %     filename = strcat(dir0(1:end-1),'.jpg');
    %     print('-djpeg','-r100',filename)
    hFig2 =figure;
    set(hFig2,'colormap',jet)
    subplot(411)
    yyaxis left
    %plot(1:numFrames,[cell_metrics.Area],'-',1:numFrames,[nuclei_metrics.Area],'--','linewidth',2)
    %ylabel('Cell(-), Nuclei(--) Areas')
%    plot(1:numFrames,[nuclei_metrics.numEndP_Skel],'-','linewidth',2)
%    ylabel('Num Endpoints Skel')
        plot(1:numFrames,[nuclei_metrics.Dist_um_s],'-','linewidth',2)
    ylabel('Speed [um/s]')

    axis tight;grid on
    
    yyaxis right
    ratioAreas = [cell_metrics.Area]./[nuclei_metrics.Area];
    %plot(1:numFrames,ratioAreas,'-.',1:numFrames,0.95*max(ratioAreas)*[cell_metrics.touchBorder],'r:','linewidth',2)
    %ylabel('Cell/Nuclei Ratio (-.)')
    plot(1:numFrames,[nuclei_metrics.numEndP_Skel],'-','linewidth',2)
    ylabel('Num Endpoints Skel')
%    plot(1:numFrames,[nuclei_metrics.numBranches_Skel],'-','linewidth',2)
%    ylabel('Num Branchpoints Skel')
    axis tight;grid on
    set(gca,'xlim',[1 numFrames])
    
    subplot(412)
    yyaxis left
    %plot(1:numFrames,[cell_metrics.NewArea],'-',1:numFrames,0.95*max([cell_metrics.NewArea])*[cell_metrics.touchBorder],'r:','linewidth',2)
    %ylabel('New Area')
    axis tight;grid on
    plot(1:numFrames,[nuclei_metrics.Dist_um_s],'-','linewidth',2)
    ylabel('Speed [um/s]')
    
    yyaxis right
    %plot(1:numFrames,[cell_metrics.Dist],'-','linewidth',2)
    %ylabel('Centroid Movement')
    plot(1:numFrames,[nuclei_metrics.PositionR],'-','linewidth',2)
    ylabel('Position Nuc')
    
    axis tight;grid on
    set(gca,'xlim',[1 numFrames])
    
    subplot(413)
    yyaxis left
%    plot(1:numFrames,[cell_metrics.MinorAxisLength]./[cell_metrics.MajorAxisLength],'-',1:numFrames,[cell_metrics.avTortuosity],'--','linewidth',2)
%    ylabel('Maj/Min Axis (-), Tortuo(-)')
%    plot(1:numFrames,[cell_metrics.MinorAxisLength]./[cell_metrics.MajorAxisLength],'-',1:numFrames,[nuclei_metrics.MinorAxisLength]./[nuclei_metrics.MajorAxisLength],'--','linewidth',2)
%    ylabel('Min/Maj Axis Cell (-), Nuc (--)')
    plot(1:numFrames,[cell_metrics.MinorAxisLength]./[cell_metrics.MajorAxisLength],'-','linewidth',2)
    ylabel('Min/Maj Axis Cell (-)')
    axis tight;grid on
    set(gca,'xlim',[1 numFrames])
    set(gca,'ylim',[0 1])
    
    yyaxis right
    plot(1:numFrames,[nuclei_metrics.PositionR],'-','linewidth',2)
    ylabel('Position Nuc')
    axis tight;grid on
    set(gca,'xlim',[1 numFrames])
    set(gca,'ylim',[-1 1])
    
     subplot(414)
     yyaxis left
     hold off
    plot(1:numFrames,[nuclei_metrics.MinorAxisLength]./[nuclei_metrics.MajorAxisLength],'-','linewidth',3)
    hold on
    plot([1 numFrames],[0.3 0.3],'r--','linewidth',0.5)
    ylabel('Min/Maj Axis Nuc (-)')
    axis tight;grid on
        set(gca,'xlim',[1 numFrames])
    set(gca,'ylim',[0.1 0.9])
    
    yyaxis right
    plot(1:numFrames,[nuclei_metrics.PositionR],'-','linewidth',2)
    ylabel('Position Nuc')
    axis tight;grid on
    set(gca,'xlim',[1 numFrames])
    set(gca,'ylim',[-1 1])

     
     
     
     
     
     %     plot(1:numFrames,[cell_metrics.Orientation],'-',1:numFrames,[nuclei_metrics.Orientation],'--','linewidth',2)
%     axis tight;grid on
%     ylabel('Cell, Nuc Orient')
%     
%     yyaxis right
%     plot(1:numFrames,[cell_metrics.Angle],'-',1:numFrames,[nuclei_metrics.Angle],'--','linewidth',1)
%     axis tight;grid on
%     ylabel('Cell, Nuc Angle')
%     set(gca,'xlim',[1 numFrames])
    
    set(gcf,'Position',[  64         278        1122         520]);
    
    
    %print('-djpeg','-r100',figOutName)

    

    
    v = VideoWriter(videoOutName, 'MPEG-4');
    %v = VideoWriter('LSM_170320_Timelapse_10_Separation.mp4', 'MPEG-4');
    
    % Remove any frames without data
    k2=1;
    for k=1:min(numFrames,size(F,2))
        %if
        if ~isempty(F(k).cdata)
            F2(k2).cdata    = F(k).cdata;
            F2(k2).colormap = F(k).colormap;
            k2=k2+1;
        end
        
    end
    
    %v = VideoWriter('LSM_170320_Timelapse_10_Separation.gif', 'GIF');
    open(v);
    writeVideo(v,F2);
    close(v);
    
    
    %     hFig2 =figure;
    %     set(hFig2,'colormap',jet)
    %     subplot(311)
    %     yyaxis left
    %     plot(1:numFrames,[cell_metrics.numEndP],'-','linewidth',2)
    %     axis tight;grid on
    %     ylabel('Num Endpoints')
    %     yyaxis right
    %     plot(1:numFrames,[cell_metrics.Dist],'--','linewidth',2)
    %     axis tight;grid on
    %     ylabel('Centroid Movement')
    %     set(gca,'xlim',[1 numFrames])
    %
    %
    %     subplot(312)
    %     yyaxis left
    %     plot(1:numFrames,[cell_metrics.NewArea]./[cell_metrics.Area],'-','linewidth',2)
    %     axis tight;grid on
    %     ylabel('New Area/Cell Area')
    %     yyaxis right
    %     plot(1:numFrames,[cell_metrics.Area],'--','linewidth',2)
    %     axis tight;grid on
    %     ylabel('Cell Area')
    %     set(gca,'xlim',[1 numFrames])
    %
    %
    %     subplot(313)
    %     yyaxis left
    %     plot(1:numFrames,[cell_metrics.MinorAxisLength]./[cell_metrics.MajorAxisLength],'-','linewidth',2)
    %     axis tight;grid on
    %     ylabel('Min Axis/Maj Axis')
    %
    %     %subplot(514)
    %     yyaxis right
    %     plot(1:numFrames,[nuclei_metrics.PositionS]/100,'--','linewidth',2)
    %     axis tight;grid on
    %     ylabel('Rel. Position Nuclei')
    %     set(gca,'xlim',[1 numFrames])
    %
    
end

