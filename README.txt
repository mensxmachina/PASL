Code for Pathway Activity Score Learning Algorithm (PASL) for dimensionality reduction of Gene Expression Data [1] 
-----------------------------------------------------------------------------
Apply  PASL and transform your data to the PASL's lower dimensional space using the steps (I) and (II):
-----------------------------------------------------------------------------
(I) Run PASL by adding your training data and geneset matrix to the apply_PASL.m script :
	Inputs:
	1. Training dataset 'X_train'  :
	    - rows:  samples
	    - columns:  features (probe sets)
    
	2. Geneset matrix 'G' :
 	    - rows: genesets
	    - columns: features (probe sets)	
	    - 'G' is a logical matrix where the rows correspond to membership to a geneset	

	3. Geneset names 'geneset_names' :
	    - 'geneset_names' is a string array which contains the geneset names. The i-th geneset name corresponds to the i-th row of 'G'
	
	4. 'a1' :
	   - Number of atoms at the inference phase
	
	5. 'a2' :
	   - Number of atoms at the discovery phase

	6. Further hyper-parameters:
	   - 't' : parameter which defines how many times the order of genesets will be recomputed (default value = 0.9) 
	   - 'lambda': Box-Cox normalization parameter (default value = 1/3)
	   - 'm' : Number of non-zeros per atom at the discovery phase (default value = 2000)
	   - 'verbose' : logical parameter in order to display the algorithm information during running (default value = 1)		

	Outputs:
	1. Dictionary 'D' :
	    - rows: atoms that correspond to genesets
	    - columns: features (probe sets)
        - 'D' contains atoms that directly correspond to 'G'. Each atom has non-zero coefficients only for the elements that belong in a corresponding row in 'G'.

	2. Scores matrix 'L' :
	    - rows: samples
	    - columns: PASL's newly constructed features.  The j-th column corresponds to the j-th selected-by-PASL geneset

	3. 'selected_genesets' :
	    - structure array which contains the information about the genesets that are chosen for each atom of the dictionary in the inference phase
        - 'selected_genesets.geneset_ids' : is a vector with the indexes (of 'G') of the selected genesets
        - 'selected_genesets.geneset_names' : is a string array with the names of the selected genesets 

    4. 'mu' : mean value of the train data ('X_train')
    
    5. 'sigma' : standard deviation of the train data ('X_train')

    Some details:
    1. Save your results ('D', 'L', 'selected_genesets', 'mu', 'sigma') in order to project your test data to the PASL's dictionary

-----------------------------------------------------------------------------
(II) Transform your test data to the PASL's latent space of the training data by running the script transform_by_PASL.m :

	Inputs:
	1. A test dataset 'X_test' :
	    - rows:  samples
	    - columns:  features (probe sets)

	2. The outputs of the 'apply_PASL.m' function :
 	    - the dictionary 'D'
	    - the 'selected_genesets' structure array
        - 'mu' : mean value of the train data
        - 'sigma' : standard deviation of the train data

	Outputs:
	1. Scores Matrix 'L_test' :
 	    - rows: samples
	    - columns: PASL's newly constructed features. The j-th column corresponds to the j-th geneset of 
            'selected_genesets.geneset_names' (and 'selected_genesets.geneset_ids') vector

-----------------------------------------------------------------------------


For the discovery phase, the spca.m function of SpaSM [2] is used.

[1] Karagiannaki, Ioulia, et al. "Pathway Activity Score Learning for Dimensionality Reduction of Gene Expression Data." International Conference on Discovery Science. Springer, Cham, 2020.

[2] Sjöstrand, Karl, et al. "Spasm: A matlab toolbox for sparse statistical modeling." Journal of Statistical Software Accepted for publication (2012).


contact info: ioulia.karagiannaki@gmail.com