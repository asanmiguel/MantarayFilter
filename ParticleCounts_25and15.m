% Description: Process sample collection images for lobe filtration
% experiments using 25 µm (red channel) and 15 µm (green channel) particles. 
% Obtain efficiency and concentration ratios for each experiment.

% By: Andrew Clark
% Written: 02-11-2021
% Updated: 05-12-2021

close all;
clear;
clc;

% Initialization %
q = 0;
% flow = [1 2 4 6 8 10 12 16]; % inlet flow rates - oblong lobe device
flow = [1 2 4 6 8 12 16 20]; % inlet flow rates - oblong lobe device
expName = 'BentLobe_R1'; % Experiment name where sample images are located - change for individual analysis
N = length(flow);
Inlet_Matrix_red = zeros(3,N);
out1_matrix_red = zeros(3,N);
out2_matrix_red = zeros(3,N);
Inlet_Matrix_g = zeros(3,N);
out1_matrix_g = zeros(3,N);
out2_matrix_g = zeros(3,N);

for ii = flow
    q = q+1;
    
    for jj = 1:3
    
    % Read images
    inlet_im = imread([num2str(ii) '-in-' num2str(jj) '.tif']);    
    out1_im = imread([num2str(ii) '-1-' num2str(jj) '.tif']); 
    out2_im = imread([num2str(ii) '-2-' num2str(jj) '.tif']); 
    
    % Red channel images (25 µm) first
    % Inlet counts
    inlet_red = inlet_im(:,:,1); % Obtain red channel of image
    inlet_red_bw = imbinarize(inlet_red,.2); % Binarize for counts
    inlet_red_bw = bwareaopen(inlet_red_bw,7); % Remove particles too small 
    inlet_red_count = length(regionprops(inlet_red_bw,'Area')); % Obtain counts
    Inlet_Matrix_red(jj,q) = inlet_red_count; % Store in matrix for later analysis
    % REPEAT STEPS BELOW %
    % Outlet counts
    out1_red = out1_im(:,:,1);
    out2_red = out2_im(:,:,1);
    out1_red_bw = imbinarize(out1_red,.2);
    out1_red_bw = bwareaopen(out1_red_bw,7);
    out2_red_bw = imbinarize(out2_red,.2);
    out2_red_bw = bwareaopen(out2_red_bw,7);
    out1_red_count = length(regionprops(out1_red_bw));
    out2_red_count = length(regionprops(out2_red_bw));
    out1_matrix_red(jj,q) = out1_red_count;
    out2_matrix_red(jj,q) = out2_red_count;    
    
    % Green channel images (15 µm) second
    % Inlet counts
    inlet_g = inlet_im(:,:,2); % green channel of image
    inlet_g_bw = imbinarize(inlet_g,.2);
    inlet_g_bw = bwareaopen(inlet_g_bw,2);% 15 um particles are smaller so smaller minimum pixel size
    inlet_g_count = length(regionprops(inlet_g_bw,'Area'));
    Inlet_Matrix_g(jj,q) = inlet_g_count;
    % Outlet counts
    out1_g = out1_im(:,:,2);
    out2_g = out2_im(:,:,2);
    out1_g_bw = imbinarize(out1_g,.2);
    out1_g_bw = bwareaopen(out1_g_bw,2);
    out2_g_bw = imbinarize(out2_g,.2);
    out2_g_bw = bwareaopen(out2_g_bw,2);
    out1_g_count = length(regionprops(out1_g_bw));
    out2_g_count = length(regionprops(out2_g_bw));
    out1_matrix_g(jj,q) = out1_g_count;
    out2_matrix_g(jj,q) = out2_g_count;   
    
    end
end

%% Obtain efficiency and concentration ratio for individual experiments

% Preallocate for speed
inlet_red_avg =  zeros(N,1);
inlet_green_avg = zeros(N,1);
out1_red_avg = zeros(N,1);
out1_green_avg = zeros(N,1);
out2_red_avg = zeros(N,1);
out2_green_avg = zeros(N,1);
red_eff = zeros(N,1);
green_eff = zeros(N,1);
red_cr = zeros(N,1);
green_cr = zeros(N,1);

% Obtain average count of each sample
for ii = 1:length(flow)
        inlet_red_avg(ii) =  mean(Inlet_Matrix_red(:,ii));
        inlet_green_avg(ii) = mean(Inlet_Matrix_g(:,ii));
        out1_red_avg(ii) = mean(out1_matrix_red(:,ii));
        out1_green_avg(ii) = mean(out1_matrix_g(:,ii));
        out2_red_avg(ii) = mean(out2_matrix_red(:,ii));
        out2_green_avg(ii) = mean(out2_matrix_g(:,ii));
end
% Efficiency and concentration ratios calculations
for kk = 1:length(inlet_red_avg)
    % Calculate efficiency
    red_eff(kk) = 1-out2_red_avg(kk)./inlet_red_avg(kk);
    green_eff(kk) = 1-out2_green_avg(kk)./inlet_green_avg(kk);
    % Calculate concentration ratio
    red_cr(kk) = out1_red_avg(kk)./inlet_red_avg(kk);
    green_cr(kk) = out1_green_avg(kk)./inlet_green_avg(kk);
end

%% Write data to csv for later analysis 
% dlmwrite('Inlet_Red_Counts',Inlet_Matrix_red);
% dlmwrite('Out1_Red_Counts',out1_matrix_red);
% dlmwrite('Out2_Red_Counts',out2_matrix_red);
% 
% dlmwrite('Inlet_Green_Counts',Inlet_Matrix_g);
% dlmwrite('Out1_Green_Counts',out1_matrix_g);
% dlmwrite('Out2_Green_Counts',out2_matrix_g);

%% Write efficiency data to Matlab file and CR data for plotting
% Change 1st paramter with experiment name %
save([expName,'.mat'],'red_eff','green_eff','red_cr','green_cr');

% NEXT STEPS:  Saved files are used for plotting. Plots calculate the mean efficiency at
% each inlet flow rate for N = 3 experiments. Error bars are standard
% deviation.






