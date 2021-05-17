% Purpose: Write a function for the bin range portion of the particle range
%          code. Help with visualization of particle size discretion.
% Written: 03-11-2021
% By: Andrew S Clark
% Updated: 05-12-2021

function [effAvg,effStd,binRange] = binRangeSize(in,o1,o2,lowSize,binSize,hiSize,N,flow,~)

binRange = lowSize:binSize:hiSize;
flow(1) = flow(1)-.05; % Adjust for error bar 
flow(end) = flow(end)+.05; % Adjust for error bar

% Preallocate counts
histCount_in = cell(1,N);
histCount_1 = cell(1,N);
histCount_2 = cell(1,N);

% Preallocate stats
eff = zeros(numel(flow),numel(binRange)-1);
cr = zeros(numel(flow),numel(binRange)-1);
effCell = cell(1,N);
crCell = cell(1,N);
effMat  = zeros(numel(flow),numel(binRange)-1,N);
crMat  = zeros(numel(flow),numel(binRange)-1,N);

% Loop analysis over all 3 runs
for kk = 1:N
    % Place particle sizes into bins
    for ii = 1:length(flow)
        histCount_in{kk}{ii} = histcounts(in{kk}{ii},[binRange Inf]); % Put into histogram columns
        histCount_1{kk}{ii} = histcounts(o1{kk}{ii},[binRange Inf]);
        histCount_2{kk}{ii} = histcounts(o2{kk}{ii},[binRange Inf])./2; % Samples were concentrated in half before taking images to get more particle counts
    end
        
    % Determine efficiencies and concentration ratios at each size range
    for jj = 1:length(flow)
        temp_in = histCount_in{kk}{jj}; % temp counts for individual flow
        temp_1 = histCount_1{kk}{jj};
        temp_2 = histCount_2{kk}{jj};
        if jj == 1 % all flows but 2mL were concentrated so adjust for that (see above)
            temp_2 = temp_2.*2;
        end
%         
        for tt = 1:length(binRange)-1
            if temp_1(tt)+temp_2(tt) == 0
                eff(jj,tt) = NaN; % efficiency cannot be calculated if no particles of this size
                cr(jj,tt) = NaN;
            else
                eff(jj,tt) = 1-temp_2(tt)./temp_in(tt); % calculate efficiency
                cr(jj,tt) = temp_1(tt)/temp_in(tt); % calculate concentration ratio
            end        
        end 
    end
    
    % Place stats into cell array
    effCell{kk} = eff;
    crCell{kk} = cr;
    % Place into 3D matrix for later process
    effMat(:,:,kk) = eff;
    crMat(:,:,kk) = cr;  
end

% Obtain mean and standard deviation of both eff and CR
effAvg = mean(effMat,3).*100;
effStd = std(effMat,1,3).*100;

end

