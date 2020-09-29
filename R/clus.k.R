usethis::use_package('factoextra')
#' Calculate number of clusters which cluster::pam needs
#'
#' @param df a data frame.
#'
#' @author Chand Zhang
#'
#' @examples
#' # example 1
#' df <- iris[,-5]; clus.k(df)
#'
#' # example 2
#' data(wine, package = 'rattle'); df <- wine[,-1]; clus.k(df)
#'
#' @export
clus.k <- function(df) {
  # 1 Normalize data ————————————————————————————————————————————————————————————————————
  df <- as.matrix(df); df.norm <- scale(df, center = T, scale = T)
  # E Normalize data ————————————————————————————————————————————————————————————————————


  # 2 Output final result ———————————————————————————————————————————————————————————————
  cat("Choose a k-value by average silhouette width.")
  factoextra::fviz_nbclust(df.norm, cluster::pam, method = "silhouette")
  # E Output final result ———————————————————————————————————————————————————————————————
}
