% main analysis
src_folder = fullfile(pwd,'src');
addpath(src_folder)

%% Analysis for Figure 2===================================================

figure2_framework

%==========================================================================


%% Analysis for Figure 3===================================================
% Accuracy testing for different method of classifying inside and outside
% points for a noisy spherical surface whith uniform random points inside a
% ball of varying radius

methods = {'insideOutside_raw','insideOutside'};

for i = 1:numel(methods)
    
    r_ball = figure3_accuracyTesting(methods{i},i);
    
end

%Write All accuracy statistics to text files-----------------------------

figure3_algorithm_wrtingStats

%% Analyse example data using InsideOutside with improved accuracy---------

figure3_algorithm_exampleData

%==========================================================================


%% Analysis for Figure 4===================================================

figure4_application

%==========================================================================

%% Analysis for Figure 4===================================================
% Accuracy testing for different method of classifying inside and outside
% points for a noisy spherical surface whith uniform random points inside a
% ball of varying radius

methods = {'insideOutside','convexHullMethod'};

for i = 1:numel(methods)
    
    figure5_accuracyTesting(methods{i},i);
    
end
%%
figure5_plotting

%==========================================================================
