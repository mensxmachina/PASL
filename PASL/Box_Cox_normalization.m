function [var_reduced] = Box_Cox_normalization(lambda, var_reduced, nnz_geneset)
    if lambda > 0
        var_reduced = (lambda*var_reduced)./ (nnz_geneset^lambda - 1);
    elseif lambda == 0
        var_reduced = var_reduced ./ log2(nnz_geneset);
    end
end
