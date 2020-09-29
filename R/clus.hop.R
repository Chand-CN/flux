usethis::use_package('factoextra')
#' Calculate Hopkins statistic
#'
#' @param df a data frame.
#' @param rept number of the repeat times.
#'
#' @author Chand Zhang
#'
#' @examples
#' # example 1
#' df <- iris[,-5]; clus.hop(df, 10)
#'
#' # example 2
#' data(wine, package = 'rattle'); df <- wine[,-1]; clus.hop(df, 10)
#'
#' @export
clus.hop <- function(df, rept = 30) {
  # 1 Normalize data ————————————————————————————————————————————————————————————————————
  df <- as.matrix(df); df.norm <- scale(df, center = T, scale = T)
  # E Normalize data ————————————————————————————————————————————————————————————————————


  # 2 CORE CODE: Calculate Hopkins statistic ————————————————————————————————————————————
  hopkins <- function() {
    hop <- vector(); ifelse(rept == 30, 30, rept)
    for (i in 1:rept) {
      stok <- factoextra::get_clust_tendency(df.norm, n = nrow(df)-1 , seed = i, graph = F)
      hop <- append(hop, stok$hopkins_stat)
    }
    out <- round(summary(hop), 3)
    assign('out', out, envir = parent.env(environment()))
  }
  # E CORE CODE: Calculate Hopkins statistic ————————————————————————————————————————————


  # 3 Output final result ———————————————————————————————————————————————————————————————
  hopkins()
  cat('There are Hopkins statistic for', rept, 'times:', '\n', '\n')
  print(out)
  # 3 Output final result ———————————————————————————————————————————————————————————————
}





