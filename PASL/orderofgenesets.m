function [sorted_id_G, sorted_var_G] = orderofgenesets(X, G, maxatomsperpathway, lambda)
    g = size(G, 1);  % number of genesets
    
    % Store the explained variance for each atom
    var_G = zeros(maxatomsperpathway, g);

    id_G = ones(maxatomsperpathway, g);
    id_G = cumsum(id_G, 2);
    
    for j=1:g
        % Run PCA reduced to the genes of the j-th geneset
        X_reduced = X(:, find(G(j, :)));
        [~,~,var_reduced,~,~] = pca(X_reduced, 'Centered', false);
        
        % Box-Cox normalization in the reduced PCA's variance    
        nnz_geneset = nnz(G(j, :));
        [var_reduced] = Box_Cox_normalization(lambda, var_reduced, nnz_geneset);

        if size(var_reduced, 1) < maxatomsperpathway
            var_G(1:size(var_reduced, 1),j) = var_reduced(1:size(var_reduced, 1));
        else
            var_G(:, j) = var_reduced(1:maxatomsperpathway);
        end
    end
 
    % Now sort the var_G to figure out the order of genesets
    % Turn everything to a vector
    var_G = var_G(:);
    id_G  = id_G(:);
   
    [sorted_var_G, ordering] = sort(var_G, 'descend');
    
    % The geneset indexes should have the indexes of the genesets in the
    % order they should be tried
    sorted_id_G = id_G(ordering);
end