% Purpose: Write code to analyze velocity profiles for all flows in the
% Bent Lobe device. Code obtains velocity profile and finds the closest
% local velocity maxima and the inflection point at the love with the
% highest transverse velocity.
% Written: 03-15-2021
% By: Andrew S Clark

clc;
clear;
close all;
% Get filenames in the directory (simulation data for each inlet flow rate)
filelist = dir('*');
filenames = {filelist.name};
mask = ismember(filenames,{'.','..','Figures','.DS_Store'});
filenames(mask) = []; % Remove unwanted folders in directory
exit_idx = strfind(filenames,'mL'); % Gets flow rates from file names
N = numel(filenames); % Number of files/inlet flow rates to loop over

for fileName = 1:N
    close all
    % Get file name (or flow rate) for current file 
    currName = filenames{fileName};
    flowName{fileName} = currName(1:exit_idx{fileName}-1);
    flowNumber(fileName) = str2double(flowName{fileName});

    % Start by importing csv file
    headerLines = 6;
    dataOpen =  fopen(currName);
    deLim = ',';
    solData = textscan(dataOpen,'%d %f %f %f %f %f %f %f','delimiter',deLim,'headerlines',headerLines);

    % Turn coordinates and velocities into matricies and transform coordinates into µm 
    x_coord = cell2mat(solData(:,2)).*10^6; 
    y_coord = cell2mat(solData(:,3)).*10^6-1273; % Adjusts to make bottom of main channel y = 0
    z_coord = cell2mat(solData(:,4)).*10^6;
    x_vel = cell2mat(solData(:,6));
    y_vel = cell2mat(solData(:,7));
    z_vel = cell2mat(solData(:,8));

    % Import x coordinates for lobes
    x_lobe1 = [2404.7:475:9529.7]';
    x_lobe2 = [2642.9:475:9767.9]';
    x_lobe = [x_lobe1 x_lobe2]';
    x_lobe = x_lobe(:);
    x_lobe = [x_lobe; 10004.7]; % Immediately downstream of the lobe

    xend1 = [2455.3:475:9580.3];
    xend2 = [2693.4:475:9818.4];
    xend = [xend1; xend2];
    xend = xend(:);
    xend = [xend; 10055.3]; % Immediately upstream on next lobe

    % Acceptable range to find coordinate values
    x_rng = 4;
    z_rng = 10;
    
    % Import y ranges (for main channel) for analyses
    y_rng = [202 -2];

    %% Find lobe number with largest out velocity
    yOutMax = zeros(length(x_lobe),1); % Preallocation
    % Loop over every lobe
    for ii = 1:length(x_lobe)
        yOutMaxIdx = find(y_coord >= y_rng(2)-5 & y_coord <= y_rng(2)+5 & ...
         x_coord >= x_lobe(ii)-20 & x_coord <= xend(ii)+20 & z_coord >= 30 & z_coord <= 35);
        yOutMax(ii) = min(y_vel(yOutMaxIdx)); % Finds index with 'highest' negative outward y velocity
    end
    [yMax, yMaxIdx] = min(yOutMax);
    lobeNumber = yMaxIdx; % Index for later velocity profile analysis
   
    %% Start analysis
    % Preallocation
    lobeBoundaryLayer = zeros(1,length(x_lobe));

    % finds location of values within lobe range and corresponding velocities
    lobe_idx = find(x_coord > x_lobe(lobeNumber)-x_rng & x_coord < x_lobe(lobeNumber)+x_rng & z_coord >= 26 & z_coord <= 35); 
    x_lobe_coord = x_coord(lobe_idx); %get x coordinates of vel data of interest at lobe (including outer channels)
    y_lobe_coord = y_coord(lobe_idx); % y coordinates of velocity data of interest
    x_lobe_vel = x_vel(lobe_idx); % get actual velocity values
    y_lobe_vel = y_vel(lobe_idx);
    
%     % Find statistics about velocity profiles
%         max_vel = max(x_lobe_vel);
%         max_vel_idx = find(x_lobe_vel == max_vel);
%         y_max = y_lobe_coord(max_vel_idx);
% 
% 
%     % Plot velocity profile
%         figure(1)
%         scatter(x_lobe_vel,y_lobe_coord)
%         title('X Velocity Profile at 1st Lobe','fontsize',16)
%         xlabel('Velocity (m/s)','fontsize',16)
%         ylabel('Y Coordinate (m)','fontsize',16)
%         hold on
%         plot([min(x_lobe_vel) max(x_lobe_vel)], [y_max y_max], '--k')
%         hold off
    
    %% Narrow down to main channel

    % Narrow data down to important zone in main channel
    ymax_idx = find(y_lobe_coord <= y_rng(1),1,'last');
    ymin_idx = find(y_lobe_coord >= y_rng(2),1,'first');
    ymax_rng = find(y_lobe_coord <= y_rng(1) & y_lobe_coord >= y_rng(2));
    
    % y coordinate and velocities
    y_chnl_rng = y_lobe_coord(ymax_rng);
    x_chnl_vel = x_lobe_vel(ymax_rng);
    y_chnl_vel = y_lobe_vel(ymax_rng);
    
    % Velocity profile in main channel - Plot
%     figure(2)
%     scatter(x_chnl_vel,y_chnl_rng)
%     title(sprintf('Inlet Velocity Profile at Lobe %d', lobeNumber),'fontsize',40,'fontname','Times New Roman') 
%     ax = gca;
%     ax.FontSize = 35; 
%     xlabel('Velocity (m/s)','fontsize',32,'fontname','Times New Roman')
%     ylabel('Y Coordinate (µm)','fontsize',32,'fontname','Times New Roman')
%     hold on
%     ylim([-2 202])
    
    %% Now interpolate between data points and find inflection point
    interpPts = [-2:0.5:202];
    smooth_xVel = interp1(y_chnl_rng,x_chnl_vel,interpPts,'spline');
    plot(smooth_xVel,interpPts,'-b')
    
    % Use gradient data to find inflection point
    smooth_dy = gradient(smooth_xVel,interpPts);
    smooth_d2y = gradient(smooth_dy,interpPts);
    
    % Find inflection points
    indices = find([0 diff(sign(smooth_d2y))]~=0); % locations where second derivative crosses 0 (inflection point)
    
    % Adjust if there is 1 or 0 inflection points in data (there is no
    % inflection point)
    if length(indices) < 2 || isempty(indices) == 1
        indices = [1 length(smooth_d2y)];
    end
    
    % Adjust to make sure the inflection points are the correct way
    if sign(smooth_d2y(indices(1))) < 0
        indices(1) = [];
    end
    if sign(smooth_d2y(indices(end))) > 0
        indices(end) = [];
    end
    
    % Translate index to y coordinate location
    inflectPts = interpPts(indices);
    
    % Puts empty inflection point to center of the channel
    if isempty(inflectPts) == 1
        inflectPts = [100 100];
    end
    
    % Gets distances in microns from location (adjusts for 2 um buffer)
    topPt = 202-inflectPts(end);
    botPt = inflectPts(1)+2;
    inflect_avg = mean([topPt,botPt]); % Gets mean for slight differences from interp1 function
   
    
    %% Calculate max distance locations
    lcl_max = islocalmax(smooth_xVel,'maxnumextrema',3); % should only have 3 local mal
    lcl_max = find(lcl_max);
    max_pts = interpPts(lcl_max); % Translate locations to y coordinate
    
    if length(max_pts) == 3
        topDist = 202-max_pts(3);
        botDist = max_pts(1);

    elseif length(max_pts) == 2 % if doesn't notice middle value
        max_pts = [max_pts(1) 0 max_pts(2)];     
    else % No inflection point
        max_pts = [max_pts(1) 0 202-max_pts(1)];
    end
    % Adjust for occasional asymmetry or "weak" inflection points
    maxVelPt = round(mean([max_pts(1)+2 202-max_pts(3)])); % average out for slight differences between top/bottom channel
    %% Double check inflection point with second derivative
%     %Plot 2nd derivative
%     figure(3)
%     plot(interpPts,smooth_d2y,'-k')
%     hold on
% %     plot(interpPts,smooth_dy,'-g')
%     plot([0 200],[0 0],'-r')
%     if inflect_avg == 0 %if no inflection point - make it 100 (center channel) for visualization
%         inflect_avg = 100;
%     end
%     legend('d^2v/dy','Zero Line')  

    %% Print out data point
    % Adjust for occasional asymmetry or "weak" inflection points (usually from slow flow)
    if inflect_avg > 100
        inflect_avg = 100;
    end
    if maxVelPt > inflect_avg
        maxVelPt = round(inflect_avg);
    end
    if maxVelPt > 80 % Remove points that are actually bulk velocity 
        maxVelPt = 100;
    end

    fprintf('Inflection Point: %.2f\n',round(inflect_avg));
    fprintf('Max Vel Point: %.2f\n',maxVelPt);

    % Collect data for overall plotting
    maxVelPts(fileName) = maxVelPt;
    inflectPtAll(fileName) = round(inflect_avg);
end
%% Sort to put in order with flow rate
[flowSort, flowSortIdx] = sort(flowNumber);
maxVelPtsSort = maxVelPts(flowSortIdx);
inflectSort = inflectPtAll(flowSortIdx);

% Smooth to help later plotting
flowSmooth = flowSort(1):0.25:flowSort(end);
maxVelSmooth = interp1(flowSort,maxVelPtsSort,flowSmooth,'pchip');
inflectSmooth = interp1(flowSort,inflectSort,flowSmooth,'pchip');
topSmooth = ones(1,numel(flowSmooth))*100;

% NEXT STEPS: Plot data and overlap with particle range efficiency data


    
    