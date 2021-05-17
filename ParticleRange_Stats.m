% Written: 03-09-2021
% Description: Loads .mat files that contain particle sizes at each outlet
%       and uses those counts for various statistical analysis and figure
%       production - individual experiment results
% By: Andrew S Clark
% Updated: 05-12-2021

clear;
clc
close all
% Load data
parent = pwd;
fileNames = dir('*.mat');
expName = 'OblongLobe_eff_R1';
% Loop over each individual experiment count
for ii = 1:length(fileNames)
    % Load cell arrays
    file = fileNames(ii).name;
    inlet = load(file,'inlet_circ');
    out1 = load(file,'out1_circ');
    out2 = load(file,'out2_circ');

    % Change names of loaded variables in structure to cell array
    in{ii} = inlet.inlet_circ;
    o1{ii} = out1.out1_circ;
    o2{ii} = out2.out2_circ;
end
% Loads correct flow rate for either Oblong/Bent lobe device
if numel(in{1}) == 5
    flow = [2 4 6 10 14];
    lobeName = 'Oblong Lobe';
else
    flow = [2 4 6 10 14 18];
    lobeName = 'Bent Lobe';
end

% Obtains efficiency and std for desired bin range size - a size of 5 µm was used
% for this analysis

% [effAvg,effSTD,binRange] = binRangeSize(in,o1,o2,low size,binSize,hi Size,N,flow,lobeName);
[effAvg,effStd,binRange] = binRangeSize(inlet_circ,o1,o2,10,5,30,N,flow,lobeName);

% Save for later plotting for all N = 3 experiments
save([expName,'.mat'],'binRange','effAvg','effStd','flow');

% Efficiencies were then plotted as mean +/- standard deviation 
% Data was used to compare Range and Simulation data