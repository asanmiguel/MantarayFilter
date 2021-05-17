% Description: Code to analyze particle range experiments using images of
% droplets taken by episcope. Gets particle counts and sizes for individual
% experiments.

% Code date: 01-18-2021
% By: Andrew S Clark
% Updated: 05-12-2021

%% Obtain data
close all;
clear;
clc;

% Camera setting (Hemamatsu at 1X = 13mmx13mm FOV - 10X = 1.3mmx1.3mm)
pxl = 2048;
camFOV = 1300; % microns
um_pxl = camFOV./pxl;

% Start of new written code %

% Preallocation for speed
expName = 'OblongLobe_R1';
q = 0;
flow = [2 4 6 10 14]; % Oblong Lobe inlet flow rates
% flow = [2 4 6 10 14 18]; % Bent lobe flow rates
N = length(flow);
Inlet_Matrix = zeros(3,N);
out1_matrix = zeros(3,N);
out2_matrix = zeros(3,N);

% for loop to go over each inlet flow rate and sample image
for ii = flow
    q = q+1;
    
    % loop over each sample image of N=3
    for jj = 1:3
    % Inlet analysis first
    inlet_im = imread([num2str(ii) '-in-' num2str(jj) '_MMStack_Pos0.ome.tif']); % Read in image 
    in_bw = imbinarize(inlet_im,0.1); % Binarize image 
    [c_in, r_in] = imfindcircles(in_bw,[6, 70]); % Find circles within 6 to 70 pixels in radius
    % Outlet analysis 
    out1_im = imread([num2str(ii) '-1-' num2str(jj) '_MMStack_Pos0.ome.tif']); % Read image
    out2_im = imread([num2str(ii) '-2-' num2str(jj) '_MMStack_Pos0.ome.tif']); 
    out1_bw = imbinarize(out1_im,.1); % Binarize image
    out2_bw = imbinarize(out2_im,.1);
    
    [out1_c, out1_r] = imfindcircles(out1_bw,[6 70]); % Obtain circles within predetermined pixel size range
    [out2_c, out2_r] = imfindcircles(out2_bw,[6 70]);
    
    % Place data in cell array
    in_circ{jj} = r_in;
    o1_circ{jj} = out1_r;
    o2_circ{jj} = out2_r;
   
    end
    
    % Concatonate counts and sizes from each sample image - convert to
    % diameter using pixel to µm size
    inlet_circ{q} = cat(1,in_circ{:}).*um_pxl.*2;
    out1_circ{q} = cat(1,o1_circ{:}).*um_pxl.*2;
    out2_circ{q} = cat(1,o2_circ{:}).*um_pxl.*2;
    
end

% Export counts for later analysis in ParticleRange_Stats.m
save([expName, '.mat'],'inlet_circ','out1_circ','out2_circ');

%% Statistical portion
% Loads correct flow rate for either Oblong/Bent lobe device
if numel(in{1}) == 5
    flow = [2 4 6 10 14];
    lobeName = 'Oblong Lobe';
else
    flow = [2 4 6 10 14 18];
    lobeName = 'Bent Lobe';
end

% Calls binRange function which completes analysis and returns graph of
% efficiency by inlet flow rate for 5 um bins of particles. Also returns
% efficiency average and standard deviation for that experiment.
% [effAvg,effSTD,binRange] = binRangeSize(in,o1,o2,low size,binSize,hi Size,N,flow,lobeName);
[effAvg,effStd,binRange] = binRangeSize(inlet_circ,o1,o2,10,5,30,N,flow,lobeName);









