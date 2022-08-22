#' @title check_input
#'
#' @description Checks input data. Input data should be a 3 dimensional array
#' with samples x features x replicates or
#' list of replicate matrices or dataframes. Each matrix or dataframe
#' corresponds to a replicate and
#' should have matching dimensions (samples x features).
#' Column names should correspond to feature names
#' and must match between replicates.
#'
#' @param data input data provided to FeatSeek
#' @param replicates if data is a 2 dimensional array with samples x features
#' a dataframe with 2 columns named sample and replicates indicating which
#' sample corresponds to which replicate must be provided
#'
#'
#' @return input data reshaped to 3 dimensional array with
#' samples x features x replicates
#'
#' @keywords internal
check_input <- function(data, replicates=NULL){
    # check if list was provided as input
    if (inherits(data,"list") & length(data) > 0){

        # check if all elements in list are either df or matrix
        ismat <- all(vapply(data, function(x) is.matrix(x), logical(1)))
        isdf <- all(vapply(data, function(x) inherits(x,"data.frame"),
                           logical(1)))
        if (!(ismat | isdf)) stop("Not all elements in input list are either dataframes or matrices!")

        # get colnames and check if all replicates are matching
        cnames <- lapply(data, function(x) colnames(x))
        if (all(vapply(cnames, function(x) is.null(x), logical(1)))) stop(
            "No feature names given or features not in correct dimension of data array!")
        if (length(unique(cnames))!=1) stop("Feature names are not the same for all replicates!")

        # convert data to matrix and concat to 3 dim array
        data <- lapply(data, function(x) array(as.matrix(x), dim=dim(x)))
        data <- Reduce(function(x,y) abind::abind(x,y,along=3), data)

        # set feature names of array
        dimnames(data)[[2]] <- cnames[[1]]

        if( length(dim(data)) !=3) stop("At least 2 replicates required!")

        # check if input is 3 dim array
    } else if (is.array(data)) {
        if ( length(dim(data)) ==2){
            if (is.null(replicates)) stop("Replicates not found. Please provide a dataframe indicating which sample corresponds to which replicate!")
            if (!nrow(replicates) == dim(data)[1]) stop("Replicate indicator dataframe not same length as samples in data!")
            if (!inherits(replicates, "data.frame")) stop("Please provide a dataframe indicating which sample belongs to which replicate")
            if ( ncol(replicates) != 2) stop("Please provide a dataframe indicating which sample belongs to which replicate, with columns named sample and replicate!")
            if ("sample" %in% replicates & "replicate" %in% replicates) stop("Please provide a dataframe indicating which sample belongs to which replicate, with columns named sample and replicate!")
            if (any(table(replicates$sample) < 2)) stop("Not every sample has at least 2 replicates!")
            if ( min(table(replicates$replicate)) < ncol(data)) stop("Too many missing samples for one of the replicates!")
            if( !all(replicates$sample == floor(replicates$sample)) & !all(replicates$replicate == floor(replicates$replicate))) stop("Please provide integer ids for samples and replicates!")
            nsamples <- length(unique(replicates$sample))
            nrep <-  length(unique(replicates$replicate))
            nfeat <- ncol(data)
            outdata <- array(dim=c(nsamples,nfeat, nrep))
            samplenames <- dimnames(data)[[1]]
            featnames <- dimnames(data)[[1]]

            for (i in unique(replicates$replicate)){
                outdata[ replicates$sample[replicates$replicate == i],,i] <- data[replicates$replicate == i,]
            }
            data <- outdata
            dimnames(data)<- list(samplenames, featnames,NULL)
        }else{
            # check if there are at least 2 replicates
            if( length(dim(data)) !=3) stop("At least 2 replicates required!")
            if( dim(data)[3] < 2) stop("At least 2 replicates required!")
            if (is.null(dimnames(data)[[2]]))  stop("No feature names given or features not in correct dimension of data array!")
        }

    } else{
        stop("Please provide a list of matrices/dataframes or a 3 dimensional array as input!")
    }
    data
}

#' @title init_selected
#'
#' @description Checks if preselected init features are in input data.
#' If init is NULL, it is set to feature with highest replicate correlation.
#'
#' @param init preselected starting set of features
#' @param data input data
#' (3 dimensional array with samples x features x replicates)
#'
#'
#' @return names of initial set of feature
#'
#' @keywords internal
init_selected <- function(init, data){

    features <- dimnames(data)[[2]]
    if (is.null(features)) stop("No feature names given or features not in correct dimension of data array!")
    # check if init features are in data
    if(!is.null(init) & !all(init %in% features)){
        stop("Could not find init features in data!")
    }

    # start with feature with highest replicate correlation if no init
    # features were given
    if (is.null(init)){
        cor <- apply(data, 2, function(x) {
            mean(stats::cor(x, use = "pairwise.complete.obs")
                 [upper.tri(!diag(nrow = dim(data)[3]))])
        })
        init <- names(which.max(cor))
    }
    init
}


#' @title filter
#'
#' @description This function removes features with mean pairwise replicate
#' pearson correlation < filter_thr
#'
#' @param data 3 dimensional array with samples x features x replicates.
#' @param filter_thr Mean pairwise replicate pearson correlation threshold
#'
#'
#' @return the filtered data array
#'
#' @keywords internal
filter <- function(data, init, filter_thr){
    message("Filtering out features with correlation < filter_thr
            across replicates!")
    # get mean of pairwise correlation between replicates,
    # i.e. the off-diagonal values of the correlation matrix
    cor <- apply(data,2, function(x)
        mean(stats::cor(x, use="pairwise.complete.obs")
             [upper.tri(!diag(nrow=dim(data)[3]))]))
    keep <-  cor > filter_thr
    features <- dimnames(data)[[2]]
    init_keep <- features %in% init
    features <- features[init_keep | keep]
    data <- data[, features,]
    p <- dim(data)[2] # get number of remaining features
    message("Remaining features: ", p)
    data
}

