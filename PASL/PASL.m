function [D, L, selected_genesets, mu, sigma] = ...
    PASL(X, G, geneset_names, a1, a2, t, lambda, m, verbose)
%%
% Pathway Activity Score Learning Algorithm (PASL) for dimensionality reduction of
% Gene Expression Data
% 
% [D, L, selected_genesets, mu, sigma] = ...
%     PASL(X, G, geneset_names, a1, a1, a2, t, lambda, m, verbose)
% Learns a latent representation for the data based on the input geneset matrix.
% The data X could be analysed as the dot product of L and D: X = L*D

% Inputs:
% 'X' is a row-sample matrix of the data
% 'G' is a logical matrix where the rows correspond to membership to a geneset.
% 'geneset_names' is a string array which contains the geneset names. The i-th geneset 
%       name corresponds to the i-th row of 'G'
% 'a1' is the number of dictionary atoms for the inference phase
% 'a2' is the number of dictionary atoms for the discovery phase
% 't' is a parameter which defines how many times the order of genesets will be recomputed
% 'lambda' is the Box-Cox normalization parameter (normalization of PCA's variance)
% 'm' the number of non-zeros per atom of the discovery phase
% 'verbose' is a logical parameter in order to display the algorithm information during running

% Outputs:
% 'D' is the dictionary matrix D, where the rows are the atoms
% in the dictionary. D contains atoms that directly correspond to the geneset 
% matrix passed, that is an atom has non-zero coefficients only for the
% elements that belong in a corresponding row in 'G'.
% 'L' corresponds to the scores matrix
% 'selected_genesets' is a structure array which contains the information
% about the genesets that are chosen for each atom in the inference phase. A geneset may be included several
% times in the dictionary, but with different coefficients. It has two
% fields:
%       'selected_genesets.geneset_ids' : is a vector with the indexes of the genesets in 'G'
%       'selected_genesets.geneset_names' :is a string array with the names of the genesets 
% 'mu' is the mean value of the data 'X'
% 'sigma' is the standard deviation of 'X' 
%%
if nargin < 9
  verbose = 1;
end
if nargin < 8
  m = 2000;
  fprintf('Run PASL with the default m value: %d \n', m)  
end
if nargin < 7
  lambda = 1/3;
  fprintf('Run PASL with the default lambda value: %f \n', lambda)  
end
if nargin < 6
  t = 0.9; 
  fprintf('Run PASL with the default t value: %f \n', t)  
end

warning('off', 'stats:pca:ColRankDefX')

[n, p] = size(X);
[g, pprime] = size(G);
assert(p == pprime, 'The data must have the same dimension as the geneset matrix. Samples and genesets in the rows.');
assert(t >= 0,'t should be between 0 and 1');
assert(t <= 1,'t should be between 0 and 1');
assert(lambda >= 0,'lambda should be between 0 and 1');
assert(lambda <= 1,'lambda should be between 0 and 1');
D = zeros(0, p);




% That's the maximum number of atoms one can make using only one geneset
maxatomspergeneset = min([a1, n-1, p]);

% The indexes and names of the selected genesets
sel_geneset_ids = []; sel_geneset_names = [];

% Inference Phase

[X_z, mu, sigma] = zscore(X);
X = X_z;
   
i   = 1; % running geneset index 
cnt = 1; % atom counter

% First define the order of the genesets
[sorted_id_G, sorted_var_G] = orderofgenesets(X, G, maxatomspergeneset, lambda);


while cnt <= a1
    % At each iteration add the atom that contributes the most to the
    % explanation of the variance, based on the current ordering from orderofgenesets function.
    
    % The geneset to be added to the dictionary (if it is close to the expected variance)
    geneset = sorted_id_G(i);
    
    % Run PCA restricted to the genes of the geneset
    X_reduced = X(:, find(G(geneset, :)));
    [d_reduced,~,var_reduced,~,~] = pca(X_reduced, 'Centered', false, 'NumComponents', 1);
   
    % Box-Cox normalization of the variance
    nnz_geneset = nnz(G(geneset, :));
    [var_reduced] = Box_Cox_normalization(lambda, var_reduced, nnz_geneset);
        
    % Check whether the variances are close enough. If not recompute the
    % order of genesets
    if var_reduced(1) / sorted_var_G(i) <= t 
        
        if verbose
            fprintf('Recompute the order of the genesets \n');
        end

        [sorted_id_G, sorted_var_G] = orderofgenesets(X, G, maxatomspergeneset, lambda);
        
        i = 1; % Reset the geneset index counter
        geneset = sorted_id_G(i);
       

        % Run PCA restricted to the genes in the geneset
        X_reduced = X(:, find(G(geneset, :)));
        [d_reduced,~,var_reduced,~,~] = pca(X_reduced, 'Centered', false, 'NumComponents', 1);
        
        % Box-Cox normalization of the variance
        nnz_geneset = nnz(G(geneset, :));
        [var_reduced] = Box_Cox_normalization(lambda, var_reduced, nnz_geneset);

    end
        
    % Add new atom to the dictionary
    D1                         = zeros(1, p);
    D1(1, find(G(geneset, :))) = d_reduced';
    D(cnt, :)                       = D1;

    % Keep indexes and names of the final selected genesets when new atom is added
    sel_geneset_ids = [sel_geneset_ids; geneset];
    sel_geneset_names = [sel_geneset_names; string(geneset_names(sel_geneset_ids(cnt)))];


    % Remove the contribution of the new atom from the data
    X = X - (X * D(cnt, :)') * D(cnt, :);  
    
    if verbose
        fprintf('Atom %d, geneset_id = %d : %s\n', ...
        cnt, sel_geneset_ids(cnt), string(geneset_names(sel_geneset_ids(cnt))));
    end
        
    i   = i + 1;     
    cnt = cnt + 1;
end
selected_genesets.geneset_ids = sel_geneset_ids;
selected_genesets.geneset_names = sel_geneset_names;

if a2 > 0 
    addpath('PASL/SpaSM')
    % Discovery Phase
    fprintf('Discovery phase ...\n');
    
    X = zscore(X);

    delta                 = inf; 
    [B]   = spca(X, [], a2, delta, -m);
    D2= B';

    % Merge the dictionary of the inference and the discovery phase
    D = [D; D2];
end

% Project the initial data to the new latent space using least squares
% Scores, i.e., latent variable values
L = X_z * pinv(D);

end
