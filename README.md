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
library(FeatSeekR)
set.seed(111)
# simulate data with 500 conditions, 3 replicates (1500 samples) and 5 latent factors 
# generating 50 features
conditions <- 500
latent_factors <- 5
replicates <- 3

data <- FeatSeekR::simData(conditions=conditions,n_latent_factors =latent_factors,
  replicates = replicates)

# select the top 5 features
res <- FeatSeek(data, max_features=5)

# plot a heatmap of the top 5 selected features 
FeatSeekR::plotSelectedFeatures(res)
```


## References
<a id="1">[1]</a> 
Fischer, B., Sandmann, T., Horn, T., Billmann, M., Chaudhary, V., Huber, W. and Boutros, M., 2015. A map of directional genetic interactions in a metazoan cell. Elife, 4, p.e05464.
```
