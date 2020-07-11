usethis::use_package('vegan')
usethis::use_package('factoextra')
usethis::use_package('ggplot2')
usethis::use_package('ggfortify')
usethis::use_package('ggrepel')
#' Compute k-means clustering
#'
#' @param v a numeric matrix, data frame.
#' @param comm if TRUE, data should be Community Ecology data(especially species data such as example 1).
#' @param k number of clusters.
#' @param check_overlap if TRUE, labels won't be overlap.
#'
#' @examples
#' # example 1
#' data(doubs, package = 'ade4'); spe <- doubs$fish[-8, ]
#' kms(spe, comm = T); kms(spe, comm = T, k = 4)
#'
#' # example 2(1:50; 51:100; 101:150)
#' set.seed(9); df <- iris[sample(150,70), ]; df <- df[order(df$Species), 1:4]
#' kms(df, comm = T); kms(df, comm = T, k = 3)
#'
#' # example 3(1:59; 60:130; 131:178)
#' data(wine, package = 'rattle'); set.seed(6); df <- wine[sample(178,70), ]; df <- df[, -1]
#' kms(df, comm = T); kms(df, comm = T, k = 3, check_overlap = T)
#'
#' @export
kms <- function(v, comm = F, k = NA, check_overlap = F) {
  options(warn = -1)
  # 1 Compute Euclidean distance and hclust —————————————————————————————————————————————
  kms.hc <- function() {
    if (comm == T) {
      v.norm <- vegan::decostand(v, "normalize")
      v.ch <- vegan::vegdist(v.norm, "euc")
    }
    else {
      v.norm <- scale(v, center = TRUE, scale = TRUE)
      v.ch <- dist(v.norm)
    }
    hc = hclust(v.ch, method = "ward.D2")
    assign('v.norm', v.norm, envir = parent.env(environment()))
    assign('v.ch', v.ch, envir = parent.env(environment()))
    assign('hc', hc, envir = parent.env(environment()))
  }
  # 2 Compute average silhouette width ——————————————————————————————————————————————————
  kms.si <- function() {
    print("Please choose a k-value by average silhouette width.")
    factoextra::fviz_nbclust(v.norm, kmeans, method = "silhouette")
  }
  # 3 Compute startpoints ———————————————————————————————————————————————————————————————
  kms.sta <- function(variables) {
    clus = cutree(hc, k); groups <- as.factor(clus)
    v.means <- matrix(0, ncol(v), length(levels(groups)))
    row.names(v.means) <- colnames(v)
    for (i in 1:ncol(v)) {
      v.means[i, ] <- tapply(v.norm[, i], clus, mean)
    }
    startpoints <- t(v.means)
    assign('startpoints', startpoints, envir = parent.env(environment()))
  }
  # 4 Output k-means clustering result ——————————————————————————————————————————————————
  kms.res <- function(variables) {
    library(ggfortify)
    p <- autoplot(kmeans(v.norm, centers = startpoints), v, shape = 1, frame = TRUE, frame.type = 'norm')
    if (check_overlap == T) {
      p + ggrepel::geom_text_repel(aes(label=rownames(v)), alpha = 0.6,
                                   point.padding = NA, vjust = 1.2)
    }else {
      p + geom_text(aes(label=rownames(v)), alpha = 0.6, vjust = 1.2)
    }
  }
  # 5 Output final result ———————————————————————————————————————————————————————————————
  if (is.na(k)) {
    kms.hc(); kms.si()
  }else {
    kms.hc(); kms.sta(); kms.res()
  }
}
