prob_equal_n2 <- function(rho) {
  (9 + rho) / (9 + 13*rho + 2*rho^2)
}

rho_seq <- seq(from=0, to=20, by=0.1)
prob_seq <- prob_equal_n2(rho_seq)
plot(rho_seq, prob_seq, xlab="rho", ylab="P(T1 = T2)", main="n=2")
