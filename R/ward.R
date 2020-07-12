usethis::use_package('vegan')
usethis::use_package('pvclust')
usethis::use_package('cluster')
usethis::use_package('ape')
usethis::use_package('rattle')
#' Compute Ward cluster tree
#'
#' Make sure that there is no newline after ward().
#'
#' @param v a numeric matrix, data frame.
#' @param comm if TRUE, data should be Community Ecology data(especially species data such as example 1).
#' @param k number of clusters.
#' @param label_color if TRUE, labels will be colored.
#' @param label_cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default.
#' @param label_offset distance between text and graphics(used only when label_color = F).
#'
#' @examples
#' # suggested method to read data
#' v <- readxl::read_xlsx('data.xlsx')
#' v1 <- as.matrix(v[, -1]); rownames(v1) <- v$colname
#' ward(v1)
#'
#' # example 1
#' data(doubs, package = 'ade4'); spe <- doubs$fish[-8, ]
#' ward(spe, comm = T); ward(spe, comm = T, k = 4)
#'
#' # example 2(1:50; 51:100; 101:150)
#' set.seed(9); df <- iris[sample(150,45), ]; df <- df[order(df$Species), 1:4]
#' ward(df, comm = T); ward(df, comm = T, k = 3)
#'
#' # example 3(1:59; 60:130; 131:178)
#' data(wine, package = 'rattle'); set.seed(6); df <- wine[sample(178,45), ]; df <- df[, -1]
#' ward(df); ward(df, k = 3, label_color = F)
#' # meaning of standardization based on example 3
#' df$Alcohol <- 50*df$Alcohol + 100; ward(df, k = 3, label_color = F)
#'
#' @export
ward <- function(v, comm = F, k = NA, label_color = T, label_cex = 0.7, label_offset = 0) {
  # 1 Compute Euclidean distance and Ward ———————————————————————————————————————————————
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
    assign('v.norm', v.norm, envir = parent.env(environment()))
    assign('v.ch', v.ch, envir = parent.env(environment()))
    assign('hc', hc, envir = parent.env(environment()))
  }
  # 2 Compute average silhouette width ——————————————————————————————————————————————————
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
  # 3 Output a Ward cluster tree ————————————————————————————————————————————————————————
  ward.tre <- function() {
    if (label_color == F) {
      p <- factoextra::fviz_dend(hc, cex = label_cex, k = k, horiz = T, alpha = 0.1,
                                 color_labels_by_k = FALSE)
      print(p + labs(title = "Hierarchical Cluster Tree") + theme(plot.title = element_text(hjust = 0.5)))
    }else {
      mypal=c('#768FDf', '#CD5C5C', '#20B2AA', '#FF9169',
              '#84B7DF', '#E08A9A', '#C8E7C1','#F8C98D',
              '#A9A9A9', '#696969')
      clus = cutree(hc, k)
      plot(ape::as.phylo(hc), tip.color = mypal[clus], cex = label_cex,
           label.offset = label_offset, main = 'Hierarchical Cluster Tree')
    }
  }
  # 4 Check the robustness of clustering result —————————————————————————————————————————
  ward.rb <- function() {
    v.pv <- pvclust::pvclust(t(v.norm), method.hclust = "ward.D2",
                             method.dist = "euc", parallel = TRUE)
    plot(v.pv, main = "Cluster dendrogram with AU/BP values(%)",
         ylab = "Height", xlab = "Distance: euclidean")
    lines(v.pv)
  }
  # 5 Check the clustering result ———————————————————————————————————————————————————————
  ward.ck <- function() {
    clus = cutree(hc, k)
    sil <- cluster::silhouette(cutree(hc, k = k), v.ch)
    rownames(sil) <- row.names(v)
    plot(sil, main = "", cex.names = 0.7, col = 2:(k + 1), nmax = 100)
    title("Silhouette plot")
  }
  # 6 Output final result ———————————————————————————————————————————————————————————————
  if (is.na(k)) {
    par(mar = c(5, 4, 4, 2) + 0.1)
    ward.hc(); ward.si()
    return("Please choose a k-value by average silhouette width.")
  }else {
    step <- select.list(c('Output a Ward cluster tree',
                          'Check the clustering result by silhouette coefficient',
                          'Check the robustness of clustering result by AU/BP values'),
                        multiple = T, title = 'Which step(s) do you want to take?')
    if (any(step == 'Output a Ward cluster tree')) {
      par(mar = c(1, 2, 2, 1))
      ward.hc(); ward.tre()
      par(mar = c(5, 4, 4, 2) + 0.1)
    }
    if (any(step == 'Check the clustering result by silhouette coefficient')) {
      ward.hc(); ward.ck()
    }
    if (any(step == 'Check the robustness of clustering result by AU/BP values')) {
      ward.hc(); ward.rb()
      print('Red underlines mean p(AU) >= 0.95')
    }
  }
}
