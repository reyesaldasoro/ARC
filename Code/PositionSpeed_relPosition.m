
%% Read the files that have been stored in the current folder
cd ('/Users/ccr22/OneDrive - City, University of London/Acad/ARC_Grant/Results')
baseDir                             = 'Metrics_2018_03_22/metrics/';
MetricsDir                          = dir(strcat(baseDir,'*.mat'));
numMetrics                          = size(MetricsDir,1);
%%
for k=13%: numMetrics
    clear cent* min* max* rr cc relP* jet*
    % Load the current file to be displayed
    load(strcat(baseDir,MetricsDir(k).name));
    % Discard all cases with NaN, these are mainly when the cell is out of the field
    % of view, also when the centroid is a NaNs, it is recorded as 1 value, when it
    % should be 2
    nanVect                 = ~isnan([cell_metrics.Area]);
    % get the x,y centroid and the times
    cent1                   = [cell_metrics(nanVect).Centroid];
    centx                   = cent1(1:2:end);
    centy                   = cent1(2:2:end);
    numTimePoints           = numel(centx);
    timeVec                 = 1:numTimePoints;
    % find the extent of the cell movements to create a grid to be placed at zero
    minX                    = round(min(centx));
    maxX                    = round(max(centx));
    minY                    = round(min(centy));
    maxY                    = round(max(centy));
    [rr,cc]                 = meshgrid(linspace(minX,maxX,8),linspace(minY,maxY,8));
    
    % Determine the metric to display: 
    %%%%       relative position of the nucleus with respect to the cell,
    %       0-centre, -1:0 back, +0:1 front 
    relPositionNuclei       = [nuclei_metrics(nanVect).PositionR];
    relPositionNuclei(isnan(relPositionNuclei))=0;
    zz=zeros(size(cc));
    MetricLabel = 'Relative position of Nuclei';
    zLimits = [-1 +1];
    %%%%       minorAxis/majorAxis
%     relPositionNuclei       = [nuclei_metrics(nanVect).MinorAxisLength]./[nuclei_metrics(nanVect).MajorAxisLength];
%     relPositionNuclei(isnan(relPositionNuclei))=0;
%     zz=ones(size(cc))*mean(relPositionNuclei);
%     MetricLabel = 'Minor Axis/Major Axis';
%     zLimits  = [0 1];
    
    
    % convert the NaN in zeros and filter with a median to avoid big jumps
    relPositionNuclei2      = (medfilt1(relPositionNuclei,3,'omitnan','truncate'));
    
    
    
    % Grab the velocity that has been calibrated to microns per second, remove the
    % NaN as well
    velVector               = [cell_metrics(nanVect).Dist_um_s];
    velVector(isnan(velVector))=0;
    
    % calibrate a colour vector according to the velocity, but relative to a maximum
    % velocity of 0.45 um/s.
    jet2                    =jet;
    velVector2              = 1+round(63*velVector/max(velVector));
    velVector3              = 1+round(63*velVector/0.45);
    velVector3(velVector3>64)=64;
    
    % Start the display
    figure
    clf
    subplot(121)
    hold off
    % A line through all the track
    %plot3(centx,centy,relPositionNuclei2,'b',centx,centy,relPositionNuclei2,'bo',centx(1),centy(1),relPositionNuclei2(1),'rs',centx(end),centy(end),relPositionNuclei2(end),'g^');axis ij
    plot3(centx,centy,relPositionNuclei2,'b')
    hold on
    % start with a square, end with a diamon
    plot3(centx(1),centy(1),relPositionNuclei2(1),'rs',...
        centx(end),centy(end),relPositionNuclei2(end),'kd','markersize',15,'linewidth',3);
    % add markers with different sizes and colours according to the current speed,
    % either relative (min-max of the set) or absolute (0-0.45 um/s)
    for k2=1:numTimePoints
        % Colours spread to the extent of the values of **THIS** set
        %        plot3(centx(k2),centy(k2),relPositionNuclei2(k2),'bo','MarkerFaceColor',jet2(velVector2(k2),:),'markersize',5+round(100*velVector(k2)))
        % Colours spread to the extent of the values of **ALL** sets
        plot3(centx(k2),centy(k2),relPositionNuclei2(k2),'bo','MarkerFaceColor',jet2(velVector3(k2),:),'markersize',5+round(100*velVector(k2)))
    end
    % add a mesh at zero
    mesh(rr,cc,zz,'facecolor','none')
    grid on;axis tight
    axis ij
    
    title(MetricsDir(k).name(1:end-12),'interpreter','none')
    rotate3d on
    view(-15,20)
    zlabel(MetricLabel,'fontsize',14)
    
    subplot(122)
    % Same plot but the axis spread all the data, good to see where in the picture is
    % the cell located. 
    hold off
    %plot3(centx,centy,relPositionNuclei2,'b',centx,centy,relPositionNuclei2,'bo',centx(1),centy(1),relPositionNuclei2(1),'rs',centx(end),centy(end),relPositionNuclei2(end),'g^');axis ij
    plot3(centx,centy,relPositionNuclei2,'b',centx(1),centy(1),relPositionNuclei2(1),'rs',centx(end),centy(end),relPositionNuclei2(end),'g^');axis ij
    axis([1 cell_metrics(1).rows 1 cell_metrics(1).columns zLimits])
    grid on
    view(-15,20)

    % Change the axis for the colours
    % Colours spread to the extent of the values of **THIS** set
    %    caxis([0 max(velVector)])
    % Colours spread to the extent of the values of **ALL** sets
    caxis([0 max(0.45)])
    colormap(jet2)
    hColorbar = colorbar;
    %hColorbar.Position(4) = hColorbar.Position(4)-0.1 ;
    hBox =annotation(gcf,'textbox',...
        [hColorbar.Position(1)+0.02 0.066+hColorbar.Position(3)+hColorbar.Position(4) 0.17 0.07],...
        'String',{'Dist [um/s]'},'Edgecolor','none','fontsize',12);
    
    % Printing section
     set(gcf,'position',[38   395   962   403])
%     filename = strcat(MetricsDir(k).name(1:end-12),'_Vel_Pos.jpg');
%      filename = strcat(MetricsDir(k).name(1:end-12),'_minMAJ_Pos.jpg');
%      print('-djpeg','-r200',filename)
%      close gcf
end

