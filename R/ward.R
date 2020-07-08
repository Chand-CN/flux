usethis::use_package('vegan')
usethis::use_package('cluster')
usethis::use_package('ape')
usethis::use_package('rattle')
#' Compute Ward cluster tree
#'
#' @param v a numeric matrix, data frame.
#' @param comm if TRUE, data should be Community Ecology data(especially species data such as example 1).
#' @param k number of clusters.
#' @param cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default.
#' @param os distance between text and graphics.
#'
#' @examples
#' # suggested method to read data
#' v <- readxl::read_xlsx('data.xlsx')
#' v1 <- as.matrix(v[, -1]); rownames(v1) <- v$colname
#' flux::ward(v1)
#'
#' # example 1
#' data(doubs, package = 'ade4'); spe <- doubs$fish[-8, ]
#' ward(spe, comm = T); ward(spe, comm = T, k = 4)
#'
#' # example 2(1:50; 51:100; 101:150)
#' set.seed(9); df <- iris[sample(150,45), ]; df <- df[order(df$Species), 1:4]
#' ward(v = df, comm = T); ward(v = df, comm = T, k = 3)
#'
#' # example 3(1:59; 60:130; 131:178)
#' data(wine, package = 'rattle'); set.seed(6); df <- wine[sample(178,45), ]; df <- df[, -1]
#' ward(v = df); ward(v = df, k = 3)
#' # meaning of standardization based on example 3
#' df$Alcohol <- 100*df$Alcohol
#' ward(v = df); ward(v = df, k = 3)
#'
#' @export
ward <- function(v, comm = F, k = NA, cex = 0.7, os = 0) {
  # 1 Compute Euclidean distance and Ward————————————————————————————————————————————————
  ward.hc <- function() {
    if (comm == T) {
      v.norm <- vegan::decostand(v, "normalize")
      v.ch <- vegan::vegdist(v.norm, "euc")
    }
    else {
      v.norm <- scale(v, center = TRUE, scale = TRUE)
      v.ch <- dist(v.norm)
    }
    hc = hclust(v.ch, method = "ward.D2")
    assign('hc', hc, envir = parent.env(environment()))
    assign('v.ch', v.ch, envir = parent.env(environment()))
  }
  # 2 Compute average silhouette width———————————————————————————————————————————————————
  ward.si <- function() {
    Si <- numeric(nrow(v))
    for (k in 2:(nrow(v) - 1)){
      sil <- cluster::silhouette(cutree(hc, k = k), v.ch)
      Si[k] <- summary(sil)$avg.width
    }
    plot(1:nrow(v),Si,type = "h",main = "Silhouette-optimal number of clusters",
         xlab = "k (number of clusters)",ylab = "Average silhouette width")
    axis(1, which.max(Si), paste(which.max(Si)), font = 2)
  }
  # 3 Output a Ward cluster tree—————————————————————————————————————————————————————————
  ward.tre <- function() {
    mypal=c('#768FDf', '#CD5C5C', '#20B2AA', '#FF9169',
            '#84B7DF', '#E08A9A', '#C8E7C1','#F8C98D',
            '#A9A9A9', '#696969')
    clus=cutree(hc, k)
    plot(ape::as.phylo(hc), tip.color = mypal[clus], cex = cex, label.offset = os,
         main = 'Ward Hierachical Clustering Analysis')
  }
  # 4 Output final result and adjust figure margins——————————————————————————————————————
  if (is.na(k)) {
    par(mar = c(5, 4, 4, 2) + 0.1)
    ward.hc(); ward.si()
    return("Please choose a k value based on Average silhouette width")
  }else {
    par(mar = c(1, 2, 2, 1))
    ward.hc(); ward.tre()
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
}
