rescale <- function (nchar, low, high){
  min_d <- min(nchar)
  max_d <- max(nchar)
  rscl <- ((high - low) * (nchar - min_d))/(max_d - min_d) + low
}
