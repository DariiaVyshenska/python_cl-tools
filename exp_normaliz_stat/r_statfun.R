wx <- function(x, y){
  all_output <- wilcox.test(x,y)
  return(all_output$p.value)
}

wt <- function(x, y){
  all_output <- t.test(x,y)
  return(all_output$p.value)
}

fdr <- function(x){
  library(stats)
  p.adjust(x, method = "fdr")
}