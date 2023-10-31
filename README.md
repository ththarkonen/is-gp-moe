# Importance sampled mixtures of Gaussian process experts
Implementation of importance sampling mixtures of Gaussian process experts. This is a repository for the software used to create some of the results for the paper 
"Mixtures of Gaussian Process Experts with SMC2'' ([arxiv.org/abs/2208.12830](https://arxiv.org/abs/2208.12830)).

The software approximates posterior distributions for parameters of mixtures of Gaussian process experts. The approximation begins with randomized partitions of data which are modelled as independent Gaussian processes.
Maximum a posteriori (MAP) estimates of the Gaussian process likelihoods are used to approximate marginal likelihoods over the Gaussian process parameters. Finally, the particles are weighted and sampled according to the MAP estimates
which are then used to construct predictive distributions for the data.

# Installation
Clone or download the repository and add the include folder and subfolders to your MATLAB path. This should happen automatically on running the main script in the project root folder.

# References
If you find the software useful, please cite [arxiv.org/abs/2208.12830](https://arxiv.org/abs/2208.12830).
