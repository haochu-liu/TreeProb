P <- arg_transition_matrix(3, rho = 1)
P_hat <- transition_hat(3, P)
eig <- eigen(P_hat)

eig$values

V <- eig$vectors         # Columns are eigenvectors
D <- diag(eig$values)    # Diagonal matrix of eigenvalues

