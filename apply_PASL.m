clear; close all; clc;
addpath('PASL')

% Load data (row-sample matrix)
load('example_data/X_train.mat');

% Load the geneset matrix -- each row corresponds to a geneset
load('example_data/G.mat');

% Load the geneset names of the geneset matrix -- the i-th geneset name corresponds to the i-th row of G
load('example_data/geneset_names.mat');


% Initialization parameters
a1      = 350;  % Number of atoms (inference phase) -- Add value
a2      = 150;  % Number of atoms (discovery phase) -- Add value
t       = 0.9;  % threshold for reordering the genesets
lambda  = 1/3;  % Box - Cox normalization parameter
m       = 2000; % Number of non-zeros per atom of discovery phase
verbose = 1;

fprintf('Data and Genesets loaded. \n');
fprintf('Apply PASL.  \n');

[D, L, selected_genesets, mu, sigma] = ...
PASL(X_train, G, geneset_names, a1, a2, t, lambda, m, verbose);


% % save your results in order to use them to transform test data to the
% % PASL's latent space
% save('example_data/PASL_results.mat', 'D', 'L', 'selected_genesets', 'mu', 'sigma');




