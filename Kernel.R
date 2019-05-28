# Mallows kernel function ----
countTies <- function(x)
{
  # @param x Must be a matrix!
  # @return A vector representing number of tied pairs in each rows of "x".
  
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  } else if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  n1 <- apply(x, 1, function(u){
    tab <- table(u)
    sum(tab * (tab-1) / 2)
  })
  
  return(n1)
}

#bandwidth; r=x[i, ]; seqs=x[j, ]
Mallows_kernel = function (bandwidth, r, seqs = NULL) 
{
  if (is.vector(r)) {
    r <- matrix(r, nrow = 1)
  } else if (!is.matrix(r)) {
    r <- as.matrix(r)
  }
  kmat <- Kendall_kernel(x = r, y = seqs)
  if (is.null(seqs)) {seqs <- r}
  stopifnot(ncol(r) == ncol(seqs))
  n0 <- choose(ncol(r), 2)
  v1 <- sqrt(n0 - countTies(r))
  v2 <- sqrt(n0 - countTies(seqs))
  
  stopifnot(nrow(kmat) == length(v1))
  stopifnot(ncol(kmat) == length(v2))
  kmat <- sweep(sweep(kmat, 1, v1, "*"), 2, v2, "*")
  dmat <- sweep(sweep(-2 * kmat, 1, v1 * v1, "+"), 2, v2 * 
                  v2, "+")
  return(exp(-bandwidth * (0.25 * dmat)))
}


# Kendall kernel function ----
Kendall_kernel <- function(x, y = NULL)
{
  if (is.null(y) && is.vector(x)) {
    kmat <- kendall_corr(x, x) # kmat is number
  } else if (is.null(y) && !is.vector(x)) {
    x <- as.matrix(x)
    kmat <- kendall_corr(t(x)) # kmat is matrix
    dimnames(kmat) <- list(rownames(x), rownames(x))
  } else if (is.vector(y) && is.vector(x)) {
    kmat <- kendall_corr(x, y) # kmat is number
  } else {
    # Convert both x and y to matrix
    if (is.vector(x)) {
      x <- matrix(x, nrow = 1)
    } else if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    if (is.vector(y)) {
      y <- matrix(y, nrow = 1)
    } else if (!is.matrix(y)) {
      y <- as.matrix(y)
    }
    # Generate kernel matrix
    if (ncol(x) != ncol(y)) 
      stop("\"x\" and \"y\" need to have the same number of ranked items!")
    kmat <- kernmat(kf = kendall_corr, mat1 = x, mat2 = y) # kmat is matrix
  }
  return(kmat)
}

kendall_corr <- function(x, y = NULL)
{
  if (requireNamespace("pcaPP", quietly = TRUE)) {
    pcaPP::cor.fk(x = x, y = y)
  } else {
    stats::cor(x = x, y = y, use = "everything", method = "kendall")
  }
}

kernmat <- function(kf, mat1, mat2)
{
  # Alternative function to \code{\link{kernlab::kernelMatrix}} for two matrices \code{mat1,mat2} with observations in rows
  # @note This function uses double for loops in R and is less cumbersome when "\code{mat2}" is not equal to "\code{mat1}"
  
  if (requireNamespace("kernlab", quietly = TRUE)) {
    res1 <- kernlab::kernelMatrix(kf, mat1, mat2)
  } else {
    stopifnot(is.matrix(mat1) && is.matrix(mat2))
    res1 <- matrix(0, nrow = nrow(mat1), ncol = nrow(mat2))
    dimnames(res1) <- list(rownames(mat1), rownames(mat2))
    for (i in seq(nrow(mat1))) {
      for(j in seq(nrow(mat2))) {
        res1[i,j] <- kf(mat1[i, ], mat2[j, ])
      }
    }
  }
  return(res1)
}
