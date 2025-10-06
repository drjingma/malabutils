#' @importFrom magrittr %>%
NULL

#' Log transform positive counts
#'
#' Replaces non-zero entries with their natural log.
#' @param x Numeric vector.
#' @return Numeric vector with log-transformed non-zero entries.
#' @export
rlog <- function(x) {
  y <- x[x != 0]
  x[x != 0] <- log(y)
  return(x)
}


#' Double center a matrix using Gower centering
#'
#' @param M Numeric matrix of dimension n x n.
#' @return A double-centered matrix.
#' @export
GowerCentering <- function(M) {
  d <- dim(M); n <- d[1]
  II <- matrix(1, n, n)
  (diag(n) - (1/n) * II) %*% M %*% (diag(n) - (1/n) * II)
}

#' Geometric mean of positive values
#'
#' Computes the geometric mean of positive and non-NA entries.
#' @param x Numeric vector.
#' @param na.rm Logical; remove NA values?
#' @return Geometric mean.
#' @export
gm.mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm = na.rm) / length(x))
}


#' Geometric mean of positive components only
#'
#' @param x Numeric vector.
#' @return Geometric mean of positive entries.
#' @export
gm.mean.pos <- function(x) {
  y <- x[x > 0 & !is.na(x)]
  if (length(y) == 0) {
    stop("The input vector is all zero!")
  }
  exp(sum(log(y)) / length(y))
}


#' Centered log-ratio transformation
#'
#' @param x Numeric vector.
#' @param base Numeric; log base.
#' @param use.pseudo Logical; add 1 for zero handling?
#' @return Transformed vector.
#' @export
clr <- function(x, base = exp(1), use.pseudo = TRUE) {
  if (use.pseudo) {
    x <- x + 1
    x <- log((x / gm.mean(x)), base)
  } else {
    x <- log((x / gm.mean(x)), base)
    x[!is.finite(x) | is.na(x)] <- 0
  }
  x
}

#' Procrustes analysis between two matrices
#'
#' @param A Numeric matrix.
#' @param B Numeric matrix.
#' @return List with normalized matrices, rotation, transformed B, and RSS.
#' @export
procrustes <- function(A, B) {
  A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A.centered, type = "F") / (ncol(A) * nrow(A))
  A.normalized <- A.centered / A.size

  B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B.centered, type = "F") / (ncol(B) * nrow(B))
  B.normalized <- B.centered / B.size

  svd.res <- svd(B.normalized %*% t(A.normalized))
  U <- svd.res$u; V <- svd.res$v; RT <- V %*% t(U)

  B.transformed <- RT %*% B.normalized
  RSS <- norm(A.normalized - B.transformed, type = "F")

  list(
    A.normalized = A.normalized,
    B.normalized = B.normalized,
    rotation.mtx = RT,
    B.transformed = B.transformed,
    RSS = RSS
  )
}


#' Robust centered log-ratio transform
#'
#' @param x Numeric vector.
#' @param base Numeric; log base.
#' @return Transformed vector.
#' @details
#' This function transform nonzero counts using the centered log ratio while leaving zero counts unmodified.
#'
#' @export
rclr <- function(x, base = exp(1)) {
  gm_pos <- sum(log(x[x > 0 & !is.na(x)], base)) / sum(x > 0 & !is.na(x))
  x <- log(x, base) - gm_pos
  x[!is.finite(x) | is.na(x)] <- 0
  x
}

#' Modified centered log-ratio with epsilon
#'
#' @param mat Numeric matrix.
#' @param base Numeric; log base.
#' @param const Numeric; constant to add.
#' @return Transformed matrix.
#' @details
#' The usual CLR transformation adds a pseudocount to all zeros or all counts, causing artificial inflation of negative values. By contrast, the modified CLR transformation transforms the positive counts and ensures that the original positive counts are all greater than zero.
#'
#' @export
clr.epsilon <- function(mat, base = exp(1), const = NULL) {
  index <- which(mat > 0 & !is.na(mat))
  a <- apply(mat, 1, function(x) {
    y <- x[x > 0 & !is.na(x)]
    y <- log(y / gm.mean.pos(y), base)
    x[x > 0 & !is.na(x)] <- y
    x
  })
  a <- t(a)
  if (is.null(const)) const <- 0.01
  mat[index] <- a[index] + abs(min(a)) + const
  mat
}

#' Function to filter taxa by relative abundances and prevalence
#' @import tibble
#' @import dplyr
#' @param otu_table a matrix with rows being samples
#' @param min_abundance threshold for minimum relative abundances
#' @param min_sample threshold for minimum prevalence
#' @return A matrix of relative abundances containing the filtered features
#' @details
#' This function takes an OTU count matrix and returns a matrix of relative abundances with
#' only features that pass the threshold.
#'
#' @export
filter_by_min_abundance_min_sample <- function(otu_table, min_abundance, min_sample){
  count_table_rel <- sweep(otu_table, 1, rowSums(otu_table), FUN='/')
  relative_abundances = data.frame(count_table_rel) %>% tibble::rownames_to_column("ID")
  keep = relative_abundances %>%
    dplyr::select(-ID) %>%
    dplyr::mutate(dplyr::across(everything(), ~ . > min_abundance)) %>%
    colSums(na.rm = TRUE) > (nrow(relative_abundances) * min_sample)
  df = relative_abundances %>% dplyr::select(c(ID,names(keep)[unname(keep)]))
  df %>% tibble::column_to_rownames("ID")
}

