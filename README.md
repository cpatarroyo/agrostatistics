# agrostatistics

## General description
`agrostatistics` is an R package that contains a toolkit to infer the sexual reproduction rate of a polyploid population in a spatially explicit context. The workflow required for such inference from a real population uses different functions required to simulate the populations used for the inference and to actually infer the probability distribution of the sexual reproduction rate in a population through Approximate Bayesian Computation [^1][^2].
This approximation is applicable to polyploid populations sampled in a single timepoint and characterized by microsatellites. 

[^1]: Csilléry, K., François, O., & Blum, M. G. B. (2012). Abc: An R package for approximate Bayesian computation (ABC). Methods in Ecology and Evolution, 3(3), 475–479. https://doi.org/10.1111/j.2041-210X.2011.00179.x.
[^2]: Beaumont, M. A., Zhang, W., & Balding, D. J. (2002). Approximate Bayesian Computation in Population Genetics. Genetics, 162(4), 2025–2035. https://doi.org/10.1093/genetics/162.4.2025.
[^3]: R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

### Instalation

To install this package use: 
```
install.packages("devtools")
devtools::install_github("cpatarroyo/agrostatistics")
```

### Recomendations
As with most ABC applications, the typical pipeline for the inference of the sexual reproduction rate in a real population includes two main steps:
1. Simulation of populations and calculation of the summary statistics. 
2. Inference of the probability distribution of the parameter of interest.

The first step is the most computationally demanding process in the workflow. This is done using the `refTabMake` function that takes advantage of the `parallel` R package [^3] to use multiple cores in parallel to make the simulation step more efficient. Despite the attempts to increase efficiency, the simulation runs may take more than 24 hours depending of the amount of individuals and generations to be simulated per population. **Example** A run of 200 populations for 10.000 generations and 1.000 individuals takes around 30 hours on 10 CPUs. For this reason, it is recommended to run the simulations in a cluster or a remote server if possible. Otherwise, be advised that the computer used to run the simulations might be busy for days at a time.

**Warning**: Windows does not support multi-core execution using the `parallel` package. To use the `refTabMake` function in Windows set the `cores` parameter to 1.

On the other hand, the inference of the probability distribution of the parameter and the cross-validation take only a few minutes to execute. 
