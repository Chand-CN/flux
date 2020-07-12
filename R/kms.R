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
#' @param label_cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default.
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
#' kms(df); kms(df, k = 3)
#'
#' @export
kms <- function(v, comm = F, k = NA, label_cex = 3.5, label_alpha = 0.6, label_hjust = -0.3) {
  options(warn = -1)
  # 1 Compute Euclidean distance and hclust —————————————————————————————————————————————
  kms.hc <- function() {
    if (comm == T) {
      v.norm <- vegan::decostand(v, "normalize")
    }
    else {
      v.norm <- scale(v, center = TRUE, scale = TRUE)
    }
    assign('v.norm', v.norm, envir = parent.env(environment()))
  }
  # 2 Compute average silhouette width ——————————————————————————————————————————————————
  kms.si <- function() {
    print("Please choose a k-value by average silhouette width.")
    factoextra::fviz_nbclust(v.norm, kmeans, method = "silhouette")
  }
  # 3 Output k-means clustering result ——————————————————————————————————————————————————
  kms.res <- function() {
    library(ggfortify)
    p <- autoplot(kmeans(v.norm, centers = k, nstart = 25), v,
                  shape = 1, frame = TRUE, frame.type = 'norm')
    style <- select.list(c('points and labels - default', 'points and labels - check overlap', 'only labels'),
                        title = 'Which kind of plot do you want to choose?')
    if (style == 'points and labels - default') {
      p + geom_text(aes(label=rownames(v)), size = label_cex, alpha = label_alpha, hjust = label_hjust)
    }else if (style == 'points and labels - check overlap') {
      p + ggrepel::geom_text_repel(aes(label=rownames(v)), point.padding = NA,
                                   size = label_cex, alpha = label_alpha, hjust = label_hjust)
    }else {
      p + geom_label(label=rownames(v), size = label_cex, alpha = label_alpha)
    }
  }
  # 4 Output final result ———————————————————————————————————————————————————————————————
  if (is.na(k)) {
    kms.hc(); kms.si()
  }else {
    kms.hc(); kms.res()
  }
}
