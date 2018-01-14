bridges <- function(dat, mode = graph, connected = c("strong", "weak")) {
  e_cnt <- network::network.edgecount(dat)
  if (mode == "graph") {
    cmp_cnt <- components(dat)
    b_vec <- rep(FALSE,e_cnt)
    for(i in 1:e_cnt){
      dat2 <- dat
      delete.edges(dat2, i)
      b.vec[i] <- (components(dat2) != cmp_cnt)
    }
    
  }
  else {
    cmp_cnt <- components(dat, connected = connected)
    b_vec <- rep(FALSE, e_cnt)
    for(i in 1:e_cnt) {
      dat2 <- dat
      delete.edges(dat2, i)
      b_vec[i] <- (components(dat2, connected = connected) != cmp_cnt)
    }
  }
return(b_vec)
}
