usethis::use_package('pwr')
usethis::use_package('car')
#' Student's t-Test
#'
#' Performs one and two sample t-tests on vectors of data.  Displays results of normality test and two sample var.test.
#'
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param paired a logical indicating whether you want a paired t-test.
#' @param premise a logical indicating whether you want know results of normality test and two sample var.test are displayed.
#'
#' @examples
#' # example 1
#' set.seed(3); v1 <- rnorm(5000, 0, 1); v2 <- rnorm(5000, 0.02, 1)
#' ttest(v1 = v1, v2 = v2, alternative = 'less')
#'
#' # example 2
#' v1 <- c(2.41, 2.90, 2.75, 2.23, 3.67, 4.49, 5.16, 5.45, 2.06, 1.64, 1.06, 0.77)
#' v2 <- c(2.80, 3.04, 1.88, 3.43, 3.81, 4.00, 4.44, 5.41, 1.24, 1.83, 1.45, 0.92)
#' ttest(v1 = v1, v2 = v2, alternative = 'greater', paired = T, premise = T)
#'
#' @export
ttest <- function(v1,v2 = NA,mu = NA,alternative,paired = F,premise = F) {
  # 1 'nort' for 'normality test','tsv' for 'two sample var.test' ———————————————————————
  nort <- function() {
    ifelse(is.na(mu),rnort <- sapply(list(v1,v2), shapiro.test),rnort <- shapiro.test(v1))
    assign('rnort',rnort,envir = parent.env(environment()))}
  tsv <- function() {
    if(rnort[2]>0.05 & rnort[6]>0.05){
      rtsv <- var.test(v1,v2)
    }else{
      v<-c(v1,v2);l1<-length(v1);l2<-length(v2)
      df<-data.frame(v,F = factor(rep(1:2,c(l1,l2))))
      rtsv <- car::leveneTest(df$v, df$F)
    }
    assign('rtsv',rtsv,envir = parent.env(environment()))}
  # 2 'ost' for 'one sample t-test','tst' for 'two sample t-test'————————————————————————
  ost <- function() {t.test(v1,mu = mu,alternative = alternative)}
  tst <- function() {t.test(v1,v2,alternative = alternative,paired = paired,var.equal = T)}
  # 3 'pwrt' for 'power of t-test' ——————————————————————————————————————————————————————
  pwrt <- function() {
    l1<-length(v1);l2<-length(v2)
    d <- ifelse(is.na(mu),(mean(v1)-mean(v2))/sqrt(((l1-1)*var(v1)+(l2-1)*var(v2))/(l1+l2-2)),
                (mean(v1)-mu)/sd(v1))
    ifelse(is.na(v2),rpwrt <- pwr::pwr.t.test(n=length(v1),d=d,type='one.sample',alternative=alternative),
           ifelse(length(v1)==length(v2),ifelse(paired==F,
                                                rpwrt <- pwr::pwr.t.test(n=length(v1),d=d,type='two.sample',alternative=alternative),
                                                rpwrt <- pwr::pwr.t.test(n=length(v1),d=d,type='paired',alternative=alternative)),
                  rpwrt <- pwr.t2n.test(n1=length(v1),n2=length(v2),d=d,alternative=alternative)))
    assign('rpwrt',rpwrt,envir = parent.env(environment()))
  }
  # 4 Output final result ———————————————————————————————————————————————————————————————
  if (is.na(mu)) {
    nort();tsv()
    ifelse(rtsv[[3]]<0.05,ans <- 'Failed to pass the homogeneity test',
           ans <- list(tst(),pwrt()))
    print(ans)
  }else{
    print(list(ost(),pwrt()))
  }
  # 5 Premise(optional) —————————————————————————————————————————————————————————————————
  if (premise == T)
    if (is.na(mu)) {
      nort();tsv();list(rnort,rtsv)
    }else {nort();rnort}
}
