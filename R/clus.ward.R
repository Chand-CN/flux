usethis::use_package('factoextra')
usethis::use_package('data.table')
usethis::use_package('mclust')
usethis::use_package('ggdendro')
usethis::use_package('ggplot2')
usethis::use_package('ape')
#' Calculate Ward cluster tree
#'
#' @param df a data frame.
#' @param k number of clusters.
#' @param writetree if TRUE, output the tree to working directory.
#'
#' @author Chand Zhang
#'
#' @examples
#' # example 1
#' df <- iris[,-5]; clus.ward(df, 3)
#' table(df.gauss, iris$Species)
#'
#' # example 2
#' data(wine, package = 'rattle'); df <- wine
#' clus.ward(df[,-1], k = 3)
#' table(df.ward, df$Type)
#'
#' @export
clus.ward <- function(df, k, writetree = F) {
  # 1 Test data —————————————————————————————————————————————————————————————————————————
  df_test <- function() {
    if (sum(sapply(df, is.numeric)) + 1 < length(df)) return(cat('There exist at least two non-numeric variables.',
                                                                 '\n', 'Delete all, or keep one as cluster labels.', '\n'))
    data.table::setDT(df)
    if (sapply(df, is.character)[1] == T & length(table(df[,1])) == nrow(df)) {
      assign('del.mark', T, envir = parent.env(environment()))
    }
  }
  # E Test data —————————————————————————————————————————————————————————————————————————


  # 2 Normalize data ————————————————————————————————————————————————————————————————————
  df_norm <- function() {
    if (exists("del.mark")) {
      remove('del.mark', envir = parent.env(environment()))
      df.matri <- as.matrix(df[,-1])
      df.norm <- scale(df.matri, center = T, scale = T)
      names(df) <- LETTERS[1:length(df)]; rownames(df.norm) <- df$A
    }else {
      df.norm <- scale(df, center = T, scale = T)
    }
    assign('df.norm', df.norm, envir = parent.env(environment()))
  }
  # E Normalize data ————————————————————————————————————————————————————————————————————


  # 3 CORE CODE: Calculate Ward cluster tree ————————————————————————————————————————————
  df_ward <- function() {
    # Calculate Ward cluster tree ———————————————————————————————————————————————————————
    df.dismatri <- dist(df.norm); df.ward = hclust(df.dismatri, method = "ward.D2")
    require(ggplot2)
    p <- ggdendro::ggdendrogram(df.ward, rotate = TRUE, size = 4, theme_dendro = FALSE, color = "tomato") +
      labs(title = 'Ward Hierachical Clustering Analysis',x = '', y = '') +
      theme(plot.title = element_text(hjust = 0.5))
    plot(p)
    # Compare results with Gaussian mixture models(GMMs)—————————————————————————————————
    require(mclust); df.gauss <- Mclust(df.norm, G = k)
    df.wardcut <- cutree(df.ward, k = k)
    # Output results to global environment———————————————————————————————————————————————
    assign('df.ward', df.wardcut, envir = parent.env(environment()))
    assign('df.gauss', df.gauss$classification, envir = parent.env(environment()))
    # Write tree ————————————————————————————————————————————————————————————————————————
    if (writetree == T) {
      df.ward <- ape::as.phylo(df.ward); ape::write.tree(df.ward,file = "tree.txt")
    }
  }
  # E CORE CODE: Calculate Ward cluster tree ————————————————————————————————————————————


  # 4 Output final result ———————————————————————————————————————————————————————————————
  df_test(); df_norm(); df_ward()
  assign('df.ward', df.ward, envir = globalenv())
  assign('df.gauss', df.gauss, envir = globalenv())
  cat('There are results of ward(horizontal) compared with GMMs(vertical).', '\n', '\n')
  com <- table(df.ward, df.gauss); print(com)
  cat('\n', 'If you have your own classification, try this:', '\n', 'table(df.gauss, df$classification)')
  # E Output final result ———————————————————————————————————————————————————————————————
}






