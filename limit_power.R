limit_power <- function(P) {
  eig_r <- eigen(P)
  idx_r <- which.min(Mod(eig_r$values - 1))
  v <- Re(eig_r$vectors[, idx_r])
  
  eig_l <- eigen(t(P))
  idx_l <- which.min(Mod(eig_l$values - 1))
  w <- Re(eig_l$vectors[, idx_l])
  
  outer(v, w) / sum(w * v)
}
