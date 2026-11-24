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

## ===========================================================================
## Padded "box" numbering with an explicit entry formula
## ===========================================================================
## Coordinates: a = i + k (left lineages), b = j + k (right lineages), k.
## Every (a, b, k) with 1 <= a, b <= n and 0 <= k <= n gets a position
##     r = (a - 1) n (n + 1) + (b - 1)(n + 1) + k + 1,   r = 1, ..., n^2 (n + 1),
## i.e. lexicographic order in (a, b, k) with k fastest.  Positions with
## k > min(a, b) (i.e. i < 0 or j < 0) are padding, not real states; their
## rows and columns are all zero.  No real state moves to a padding position.
## Each move is a fixed shift in r, so the non-zeros lie on five diagonals:
##   Re  (k-1)              : c = r - 1,             weight k rho
##   C2  (k+1)              : c = r + 1,             weight (a-k)(b-k)
##   Cl1+Cl3 (a-1)          : c = r - n(n+1),        weight C(a,2) - C(k,2)
##   Cr1+Cr3 (b-1)          : c = r - (n+1),         weight C(b,2) - C(k,2)
##   C4  (a-1, b-1, k-1)    : c = r - n(n+1) - n - 2, weight C(k,2)
## all divided by D = C(a + b - k, 2) + k rho.  The MRCA (0,0,1) is position
## r = 2 and is absorbing; the start state (0,0,n) is the last position.

## (a, b, k) -> r   (vectorised)
box_index <- function(a, b, k, n) (a - 1) * n * (n + 1) + (b - 1) * (n + 1) + k + 1

## r -> (a, b, k, i, j, valid)   (vectorised)
box_state <- function(r, n) {
  a <- (r - 1) %/% (n * (n + 1)) + 1
  b <- ((r - 1) %% (n * (n + 1))) %/% (n + 1) + 1
  k <- (r - 1) %% (n + 1)
  data.frame(r = r, a = a, b = b, k = k, i = a - k, j = b - k,
             valid = k <= pmin(a, b))
}

## All n^2 (n+1) positions of the box, in order.
arg_states_box <- function(n) {
  stopifnot(n >= 2, n == round(n))
  s <- box_state(seq_len(n^2 * (n + 1)), n)
  rownames(s) <- sprintf("[%d,%d,%d]", s$a, s$b, s$k)   # [a,b,k] labels
  s
}

## Explicit formula for a single entry P[r, c]   (vectorised in r and c).
## Returns NA when r is a padding position.
P_entry <- function(r, c, n, rho) {
  s <- box_state(r, n); a <- s$a; b <- s$b; k <- s$k
  B <- n * (n + 1); N1 <- n + 1
  D <- choose(a + b - k, 2) + k * rho
  val <- (k * rho                         * (c == r - 1) +
            (a - k) * (b - k)               * (c == r + 1) +
            (choose(a, 2) - choose(k, 2))   * (c == r - B) +
            (choose(b, 2) - choose(k, 2))   * (c == r - N1) +
            choose(k, 2)                    * (c == r - B - N1 - 1)) / D
  mrca <- a == 1 & b == 1 & k == 1
  val[mrca] <- as.numeric(c == r)[mrca]
  val[!s$valid] <- NA
  val
}

## Full transition matrix in the box numbering, built from the formula.
##   n        : sample size (>= 2)
##   rho      : scaled recombination rate (R = 2 N r in the paper)
##   symbolic : TRUE -> character matrix of entries as functions of rho
##   sparse   : TRUE -> sparse matrix (package Matrix; ships with R).
##              Recommended for large n: at most 5 non-zeros per row.
arg_transition_matrix_box <- function(n, rho = 1, symbolic = FALSE, sparse = FALSE) {
  s <- arg_states_box(n)
  N <- nrow(s); B <- n * (n + 1); N1 <- n + 1
  a <- s$a; b <- s$b; k <- s$k
  live <- s$valid & !(a == 1 & b == 1 & k == 1)       # real, non-absorbing rows
  mrca <- box_index(1, 1, 1, n)
  
  ## the five moves: shift in r, constant and rho parts of the numerator
  moves <- list(
    Re  = list(shift = -1,          c0 = 0 * k,                         cr = k),
    C2  = list(shift = +1,          c0 = (a - k) * (b - k),             cr = 0 * k),
    L   = list(shift = -B,          c0 = choose(a, 2) - choose(k, 2),   cr = 0 * k),
    Rt  = list(shift = -N1,         c0 = choose(b, 2) - choose(k, 2),   cr = 0 * k),
    C4  = list(shift = -B - N1 - 1, c0 = choose(k, 2),                  cr = 0 * k)
  )
  d0 <- choose(a + b - k, 2); dr <- k                 # D = d0 + dr * rho
  
  from <- to <- integer(0); c0 <- cr <- numeric(0)
  for (m in moves) {
    use <- which(live & (m$c0 + m$cr) > 0)
    from <- c(from, use); to <- c(to, use + m$shift)
    c0 <- c(c0, m$c0[use]); cr <- c(cr, m$cr[use])
  }
  stopifnot(all(s$valid[to]))                          # never lands on padding
  lab <- rownames(s)
  
  if (symbolic) {
    P <- matrix("0", N, N, dimnames = list(lab, lab))
    P[cbind(from, to)] <- mapply(.frac_string, c0, cr, d0[from], dr[from])
    P[mrca, mrca] <- "1"
  } else {
    stopifnot(is.numeric(rho), length(rho) == 1, rho >= 0)
    x <- (c0 + cr * rho) / (d0[from] + dr[from] * rho)
    if (sparse) {
      P <- Matrix::sparseMatrix(i = c(from, mrca), j = c(to, mrca), x = c(x, 1),
                                dims = c(N, N), dimnames = list(lab, lab))
    } else {
      P <- matrix(0, N, N, dimnames = list(lab, lab))
      P[cbind(from, to)] <- x
      P[mrca, mrca] <- 1
    }
  }
  attr(P, "states") <- s
  P
}

## Convert a box-numbered matrix (numeric, symbolic or sparse) to the old
## compact matrix: drop the padding and reorder to the order of arg_states(n).
box_to_old <- function(Pbox) {
  N <- nrow(Pbox)
  n <- round(uniroot(function(x) x^2 * (x + 1) - N, c(1, N))$root)
  if (n^2 * (n + 1) != N) stop("dimension ", N, " is not n^2 (n+1) for any n")
  old <- arg_states(n)                                 # (i, j, k) in old order
  pos <- box_index(old$i + old$k, old$j + old$k, old$k, n)
  P <- Pbox[pos, pos, drop = FALSE]
  dimnames(P) <- list(rownames(old), rownames(old))
  attr(P, "states") <- old
  P
}

## ---- example ---------------------------------------------------------------
if (sys.nframe() == 0) {
  P2 <- arg_transition_matrix(2, symbolic = TRUE)             # = Eq. (2)
  print(noquote(unclass(P2)[, ]))
  P3 <- arg_transition_matrix(3, rho = 0.5)
  print(dim(P3)); print(range(rowSums(P3)))
  
  ## box numbering, n = 2 (12 x 12, padding rows/columns are zero)
  B2 <- arg_transition_matrix_box(2, symbolic = TRUE)
  print(noquote(unclass(B2)[, ]))
  ## back to the old matrix: identical to Eq. (2) above
  print(identical(unclass(box_to_old(B2))[, ], unclass(P2)[, ]))
}

## ---- example ---------------------------------------------------------------
if (sys.nframe() == 0) {
  P2 <- arg_transition_matrix(2, symbolic = TRUE) # = Eq. (2)
  print(noquote(unclass(P2)[, ]))
  P3 <- arg_transition_matrix(3, rho = 0.5)
  print(dim(P3)); print(range(rowSums(P3)))
}
