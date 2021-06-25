clear; close all; clc;

% Load the PASL results from the training phase
load('example_data/PASL_results.mat');

% Load your test data (row-sample matrix)
load('example_data/X_test.mat')


% Subtract the mean and divide by the standard deviation of the train data
% (the data that are used to construct the dictionary)
X_test = X_test - repmat(mu, size(X_test, 1), 1);
X_test = X_test ./ repmat(sigma, size(X_test, 1), 1);

% Transform your test data to the PASL's latent space
% The j-th column corresponds to the j-th geneset of 'selected_genesets.geneset_names'
L_test = X_test * pinv(D);