## Discrete-time (embedded) transition matrix of the two-locus (m = 2)
## tree-building Markov chain of Simonsen & Churchill (1997), Section 3.1,
## generalised to a sample of n >= 2 sequences.
##
## State (i, j, k): i = # single-left (type 1), j = # single-right (type 2),
##                  k = # double (type 3) sequences, subject to
##                  1 <= i + k <= n,  1 <= j + k <= n,  i, j, k >= 0.
## Number of states: n * (2 n^2 + 9 n + 1) / 6.
##
## Event rates (Table 2), all divided by  D = choose(i + j + k, 2) + k rho :
##   Re  : k rho        (i+1, j+1, k-1)
##   Cl1 : choose(i, 2) (i-1, j,   k  )
##   Cr1 : choose(j, 2) (i,   j-1, k  )
##   C2  : i j          (i-1, j-1, k+1)
##   Cl3 : i k          (i-1, j,   k  )
##   Cr3 : j k          (i,   j-1, k  )
##   C4  : choose(k, 2) (i,   j,   k-1)
## The MRCA state (0, 0, 1) is absorbing.
##
## State order: start state (0, 0, n) first and the MRCA (0, 0, 1) last.
## In between, states are sorted by the total number of ancestral lineages
## (i + k) + (j + k) (decreasing), then k (decreasing), then i (decreasing).
arg_states <- function(n) {
  stopifnot(n >= 2, n == round(n))
  g <- expand.grid(i = 0:n, j = 0:n, k = 0:n)
  is_state <- with(g, i + k >= 1 & i + k <= n & j + k >= 1 & j + k <= n)
  s <- g[is_state, ]
  mrca <- s$i == 0 & s$j == 0 & s$k == 1
  ord <- order(mrca, -(s$i + s$j + 2 * s$k), -s$k, -s$i)
  s <- s[ord, ]
  rownames(s) <- sprintf("(%d,%d,%d)", s$i, s$j, s$k)
  stopifnot(nrow(s) == n * (2 * n^2 + 9 * n + 1) / 6)
  s
}

## Transition matrix.
##   n        : sample size (>= 2)
##   rho      : scaled recombination rate (R in the paper) (numeric, >= 0)
##   symbolic : if TRUE, return a character matrix of the entries written
##              as functions of rho (rho is then ignored)
arg_transition_matrix <- function(n, rho = 1, symbolic = FALSE) {
  s <- arg_states(n)
  ns <- nrow(s)
  key <- rownames(s)
  idx <- setNames(seq_len(ns), key)
  
  if (symbolic) {
    P <- matrix("0", ns, ns, dimnames = list(key, key))
  } else {
    stopifnot(is.numeric(rho), length(rho) == 1, rho >= 0)
    P <- matrix(0, ns, ns, dimnames = list(key, key))
  }
  
  for (a in seq_len(ns)) {
    i <- s$i[a]; j <- s$j[a]; k <- s$k[a]
    
    if (i == 0 && j == 0 && k == 1) {
      P[a, a] <- if (symbolic) "1" else 1
      next
    }
    
    C <- choose(i + j + k, 2)
    ## event: c(numerator constant, numerator rho-coefficient, di, dj, dk)
    ev <- rbind(
      Re  = c(0,            k, +1, +1, -1),
      Cl1 = c(choose(i, 2), 0, -1,  0,  0),
      Cr1 = c(choose(j, 2), 0,  0, -1,  0),
      C2  = c(i * j,        0, -1, -1, +1),
      Cl3 = c(i * k,        0, -1,  0,  0),
      Cr3 = c(j * k,        0,  0, -1,  0),
      C4  = c(choose(k, 2), 0,  0,  0, -1)
    )
    ev <- ev[ev[, 1] + ev[, 2] > 0, , drop = FALSE]
    
    ## accumulate numerators per destination (several events may share one)
    num_c <- num_r <- numeric(0)
    for (e in seq_len(nrow(ev))) {
      to <- sprintf("(%d,%d,%d)", i + ev[e, 3], j + ev[e, 4], k + ev[e, 5])
      if (!to %in% key) stop("transition leaves state space: ", key[a], " -> ", to)
      # sum numerator constant
      num_c[to] <- (if (is.na(num_c[to])) 0 else num_c[to]) + ev[e, 1]
      # sum numerator rho-coefficient
      num_r[to] <- (if (is.na(num_r[to])) 0 else num_r[to]) + ev[e, 2]
    }
    
    for (to in names(num_c)) {
      b <- idx[[to]]
      if (symbolic) {
        P[a, b] <- .frac_string(num_c[[to]], num_r[[to]], C, k)
      } else {
        P[a, b] <- (num_c[[to]] + num_r[[to]] * rho) / (C + k * rho)
      }
    }
  }
  attr(P, "states") <- s
  P
}

## ---- helpers for the symbolic output -------------------------------------
.gcd <- function(a, b) if (b == 0) abs(a) else .gcd(b, a %% b)

.lin <- function(c0, cR) {
  rpart <- if (cR == 0) "" else if (cR == 1) "rho" else paste0(cR, "*rho")
  if (c0 == 0) return(rpart)
  if (cR == 0) return(as.character(c0))
  paste0(c0, "+", rpart)
}

.frac_string <- function(nc, nr, dc, dr) {
  if (nr == 0 && dr == 0) {
    g <- .gcd(nc, dc); nc <- nc / g; dc <- dc / g
    return(if (dc == 1) as.character(nc) else paste0(nc, "/", dc))
  }
  g <- Reduce(.gcd, c(nc, nr, dc, dr))
  nc <- nc / g; nr <- nr / g; dc <- dc / g; dr <- dr / g
  num <- .lin(nc, nr); den <- .lin(dc, dr)
  if (num == den) return("1")
  if (grepl("\\+", num)) num <- paste0("(", num, ")")
  paste0(num, "/(", den, ")")
}

## ---- example ---------------------------------------------------------------
if (sys.nframe() == 0) {
  P2 <- arg_transition_matrix(2, symbolic = TRUE) # = Eq. (2)
  print(noquote(unclass(P2)[, ]))
  P3 <- arg_transition_matrix(3, rho = 0.5)
  print(dim(P3)); print(range(rowSums(P3)))
}
