# ultrametricMat
The ultrametricMat is a package designed to conduct the Bayesian inference on the ultrametric matrices by leveraging the bijection map between the ultrametric matrices and the tree space. We refer to more details to  Yao et al. (2023+) Geometry-driven Bayesian Inference for Ultrametric Covariance Matrices.

# Manual
We wrap up the MCMC algorithm in a main function of `MCMC_Tr` with four arguments required for users to run the algorithm as follows:
* `Obsdf_`: the input data matrix;
* `betaSplt.tbl_`: the prior distribution generated from the beta-splitting prior. To obtain the beta-splitting prior, users can call the function `betaSplt.tbl.gen()` with two arguments of (1) `n` to list the number of leaves and (2) `beta_` to specify the hyper-parameter $\beta\in (-2,\infty]$;
* `iteNum`: the iteration number of MCMC algorithm;
* `burnIn`: the number of iterations to be discarded.

The `MCMC_Tr` function also provides other optional arguments for users to control the algorithm. We list some important arguments as follows:
* `tipBr.sd_` and `intBr.sd_`: positive real numbers to control the step size of the truncated normal proposal function when updating the edge lengths; the default is $0.5$ in the function;
* `init.Splt` and `init.Tip`: initial edge set and edge lengths for the algorithm; the default is the random edge set from the prior;
* `hypParam.brLen_`: the mean hyper-parameter for the prior on the edge lengths; the default is 1;
* `seed_`: the real number to control the random seed
  
The `MCMC_Tr` function generates the posterior samples in a named list object, and users can access the results by the name below:
* `llh`: the log-likelihood for the algorithm;
* `rtEdge`: a real vector showing the posterior samples for the lengths of root edge;
* `split`: a list object containing the posterior edge sets and the internal edge lengths;
* `tip`: a list object including the posterior samples of the leaf edge lengths;
* `init.Tr`: a list specifying the initial tree used in the algorithm.

## Other useful functions
The package also offers other necessary functions to build the bijection map between the ultrametric matrix and the tree structure. Specifically, we store the tree object in two different formats: (1) `phylo4` object from the `phylobase` package and (2) the edge set of the matrix. We list all functions that allow users to map bijectively between the ultrametric matrices and the tree objects.
* `phylo42Splt`
* `phylo42TrCov`
* `splt2Phylo4`
* `splt2TrCov`
* `TrCov2Splt`

