states_hat <- function(n) {
  stopifnot(n >= 2, n == round(n))
  s_nrow <- n * (n + 3) / 2 - 1
  s <- data.frame(
    i=rep(NA, s_nrow),
    j=rep(NA, s_nrow),
    k=rep(NA, s_nrow)
  )
  row_idx <- 1
  for (block in n:2) {
    for (i_j in 0:block) {
      s[row_idx, ] <- c(i_j, i_j, block-i_j)
      row_idx <- row_idx + 1
    }
  }
  s[row_idx, ] <- c(0, 0, 1)
  rownames(s) <- sprintf("(%d,%d,%d)", s$i, s$j, s$k)
  stopifnot(!any(is.na(s)))
  s
}

transition_hat <- function(n, P) {
  stopifnot(nrow(P) == n * (2 * n^2 + 9 * n + 1) / 6)
  keep_states <- states_hat(n)
  P[rownames(keep_states), rownames(keep_states)]
}

## ---- example ---------------------------------------------------------------
if (sys.nframe() == 0) {
  P2 <- arg_transition_matrix(2, symbolic = TRUE)
  P2_hat <- transition_hat(2, P2)
  print(noquote(unclass(P2_hat)[, ]))
  P3 <- arg_transition_matrix(3, symbolic = TRUE)
  P3_hat <- transition_hat(3, P3)
  print(noquote(unclass(P3_hat)[, ]))
}
