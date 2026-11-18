num_states <- function(n) {
  n * (2 * n^2 + 9 * n + 1) / 6
}

num_states(2)
num_states(3)
n_seq <- 2:20
num_seq <- num_states(n_seq)
plot(n_seq, num_seq, xlab="n", ylab="number of states")
