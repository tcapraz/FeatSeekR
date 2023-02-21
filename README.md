# FeatSeekR - Package for unsupervised feature selection

A fundamental step in many analyses of high-dimensional data is dimension 
reduction. Feature selection is one approach to dimension reduction whose 
strengths include interpretability, conceptual simplicity, transferability 
and modularity.
Here, we introduce the `FeatSeekR` package, which selects features based on 
the consistency of their signal across replicates and their non-redundancy.
It takes a 2 dimensional array (samples x features) of replicated measurements
and returns a dataframe storing the selected and features ranked by 
reproducibility.

# Installation


```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("FeatSeekR")
```


# How to run


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
