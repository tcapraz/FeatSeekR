# FeatSeekR - Package for unsupervised feature selection

A fundamental step in many analyses of high-dimensional data is dimension 
reduction. Feature selection is one approach to dimension reduction whose 
strengths include interpretability, conceptual simplicity, transferability 
and modularity.
Here, we introduce the `FeatSeekR` package, which selects features based on 
the consistency of their signal across replicates and their non-redundancy.
It takes a 2 dimensional array (features x samples) of replicated measurements
and returns a `SummarizedExperiment` object storing the selected and features ranked by 
reproducibility. This work was motivated by [[1]](#1) who devised a special case of
our current method to use it on microscopy data, but did not implement it as R package.


# Installation


```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("FeatSeekR")
```


# How to run

See vignette for more detailed examples.

```{r}
set.seed(111)
# simulate data with 500 samples, 3 replicates and 5 latent factors 
# generating 50 features
samples <- 500
latent_factors <- 5
replicates <- 3

sim <- FeatSeekR::simData(samples=samples,latent_factors =latent_factors,
                replicates = replicates)
data <- sim[[1]]
reps <- sim[[2]]

# select the top 5 features
res <- FeatSeek(data, replicates =reps, max_features=5)

# plot a heatmap of the top 5 selected features 
FeatSeekR::plotSelectedFeatures(data, res)
```


## References
<a id="1">[1]</a> 
Fischer, B., Sandmann, T., Horn, T., Billmann, M., Chaudhary, V., Huber, W. and Boutros, M., 2015. A map of directional genetic interactions in a metazoan cell. Elife, 4, p.e05464.
```
