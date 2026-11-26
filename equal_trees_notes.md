# Equal trees in the two-locus chain

## Notation

States are written $(a, b, k)$:

- $a = i + k$ is the number of lineages at the left locus,
- $b = j + k$ is the number of lineages at the right locus,
- $k$ is the number of double lineages, so $i = a - k$ and $j = b - k$.

The chain starts at $(n, n, n)$ and is absorbed at the MRCA $(1, 1, 1)$. From any state there are at most five moves. Every probability is divided by $D = \binom{a+b-k}{2} + k\rho$.

| Move | Events | New state | Weight |
|---|---|---|---|
| recombination | Re | $(a, b, k-1)$ | $k\rho$ |
| left + right join | C2 | $(a, b, k+1)$ | $(a-k)(b-k)$ |
| left locus coalesces | Cl1 + Cl3 | $(a-1, b, k)$ | $\binom{a}{2} - \binom{k}{2}$ |
| right locus coalesces | Cr1 + Cr3 | $(a, b-1, k)$ | $\binom{b}{2} - \binom{k}{2}$ |
| both loci coalesce | C4 | $(a-1, b-1, k-1)$ | $\binom{k}{2}$ |

For $n = 2$, states 1, 2, 3 and 9 of Table 1 in the paper are $(2,2,2)$, $(2,2,1)$, $(2,2,0)$ and $(1,1,1)$.

---

## 1. Equal heights vs. same coalescence times at every level

### Definitions

At each locus the tree has $n-1$ coalescence events. Write $t_L^{(m)}$ for the time at which the left locus drops from $m$ to $m-1$ lineages, and define $t_R^{(m)}$ in the same way for the right locus.

- **Equal heights** (Theorem 2's event) means only the last time agrees: $t_L^{(2)} = t_R^{(2)}$, the TMRCA.
- **Same coalescence times at every level** means all of them agree: $t_L^{(m)} = t_R^{(m)}$ for every $m = n, \dots, 2$. All node heights of the two trees then coincide, so total branch lengths are equal too. The topologies may still differ (see section 2).

### In terms of the chain

Only three moves change $a$ or $b$: a left coalescence lowers $a$, a right coalescence lowers $b$, and C4 lowers both at the same moment.

- **Same times at every level** means every coalescence is a C4. On the chain this is the condition **$a = b$ at all times**. The state $(1,1,0)$ is also excluded, since it can be reached only through states with $a \ne b$.
- **Equal heights** means the two loci reach one lineage in the same event. Since $a$ and $b$ each drop by at most 1 per step, that event must be the C4 move $(2,2,2) \to (1,1,1)$. On the chain this is the condition **never visit a state with $a = 1$ or $b = 1$ before $(1,1,1)$**. States with $a \ne b$ are allowed as long as $a, b \ge 2$.

The second set of allowed states contains the first, so
$$P(\text{same times at every level}) \le P(\text{equal heights}).$$

**Why they agree for $n = 2$.** Each locus has only one coalescence, so "every level" and "the last level" are the same thing. Any state with $a \ne b$ already has one locus at 1 lineage.

### Example for $n = 3$: equal heights, but different coalescence times

The path is $(3,3,3) \xrightarrow{\text{Re}} (3,3,2) \xrightarrow{\text{Cl3}} (2,3,2) \xrightarrow{\text{Cr3}} (2,2,2) \xrightarrow{\text{C4}} (1,1,1)$.

In the table below, each lineage is written as (left samples it carries | right samples it carries). The table follows the lineages:

| Step | Event | Lineages | $(a, b, k)$ |
|---|---|---|---|
| 0 | start | (1\|1), (2\|2), (3\|3) | (3,3,3) |
| 1 | Re on sample 1 | (1\|–), (–\|1), (2\|2), (3\|3) | (3,3,2) |
| 2, time $t_1$ | Cl3: (1\|–) joins (2\|2) | (12\|2), (–\|1), (3\|3) | (2,3,2) |
| 3, time $t_2$ | Cr3: (–\|1) joins (3\|3) | (12\|2), (3\|13) | (2,2,2) |
| 4, time $t_3$ | C4 | (123\|123) | (1,1,1) |

The resulting trees are:

- **Left tree:** ((1,2),3), with nodes at times $t_1$ and $t_3$.
- **Right tree:** ((1,3),2), with nodes at times $t_2$ and $t_3$.

Both trees have height $t_3$, so this path counts towards equal heights. The first coalescences happen at different times ($t_1 \ne t_2$), and the path passes through $(2,3,2)$ where $a \ne b$, so it is excluded from "same times at every level".

### Numbers

| | Equal heights | Same times at every level |
|---|---|---|
| $n = 3$, $\rho = 0.8$ | 0.372 | 0.291 |
| $n = 10$, $\rho = 1$ | 0.243 | 0.032 |
| $n = 20$, $\rho = 1$ | 0.240 | 0.0088 |
| $n = 40$, $\rho = 1$ | 0.239 | 0.0023 |

"Same times at every level" needs all $n - 1$ coalescences to be C4 events, so its probability goes to 0 as $n$ grows. "Equal heights" depends mainly on the last few coalescences, so it levels off.

---

## 2. Why $(i, j, k)$ cannot decide whether the trees are identical

**Identical trees** means the same coalescence times and the same topology, i.e. the same set of clades.

The state $(i, j, k)$ records how many lineages of each type exist, but not which samples each lineage carries. Two histories can follow exactly the same chain path, with exactly the same events, and still end with different trees.

### Example for $n = 3$

The path, in $(i,j,k)$, is
$$(0,0,3) \xrightarrow{\text{Re}} (1,1,2) \xrightarrow{\text{Re}} (2,2,1) \xrightarrow{\text{C2}} (1,1,2) \xrightarrow{\text{C2}} (0,0,3) \xrightarrow{\text{C4}} (0,0,2) \xrightarrow{\text{C4}} (0,0,1).$$

In $(a,b,k)$ this is $(3,3,3) \to (3,3,2) \to (3,3,1) \to (3,3,2) \to (3,3,3) \to (2,2,2) \to (1,1,1)$. Every coalescence is a C4 and $a = b$ throughout, so this path counts as "same coalescence times at every level".

After the two recombinations (on samples 1 and 2), the lineages are (1|–), (–|1), (2|–), (–|2) and (3|3). The two C2 events can pair the singles in two ways:

| | History A | History B |
|---|---|---|
| C2 events | (1\|–)+(–\|1), then (2\|–)+(–\|2) | (1\|–)+(–\|2), then (2\|–)+(–\|1) |
| Doubles after C2 | (1\|1), (2\|2), (3\|3) | (1\|2), (2\|1), (3\|3) |
| First C4 | (1\|1)+(3\|3), giving (13\|13) | (1\|2)+(3\|3), giving (13\|23) |
| Left tree | ((1,3),2) | ((1,3),2) |
| Right tree | ((1,3),2) | **((2,3),1)** |
| Identical trees? | yes | **no** |

Even within History B, the outcome depends on which pair of doubles merges first. If (1|2) merges with (2|1), the result is (12|12): both loci form the clade {1,2} at the same time, and the trees are identical.

### Why the chain loses this

The chain combines the cases above into one transition:

- **C2:** the weight $ij$ counts every left-single/right-single pairing together.
- **C4:** the weight $\binom{k}{2}$ counts every pair of doubles together.

Whether the trees are identical depends on which pairing occurred, and that is not a function of the $(i,j,k)$ path.

### How big the difference is

These estimates come from the same simulated runs of the full sample-labelled process, with $\rho = 0.8$ and 20,000 runs:

| | Same times at every level | Identical trees |
|---|---|---|
| $n = 3$ | 0.2863 | 0.2856 |
| $n = 4$ | 0.2013 | 0.2000 |

The gap is small but real.

To compute identical-tree probabilities exactly, the state must record each lineage as a pair (left sample set, right sample set). That state space grows exponentially in $n$, so simulation is the practical route.

---

## 3. Eigenvalues of $\hat P$

Order the allowed states as the transient states $T$ first and the MRCA $M = (1,1,1)$ last. Then
$$\hat P = \begin{pmatrix} Q & \mathbf b \\ \mathbf 0^\top & 1 \end{pmatrix}, \qquad Q = \hat P[T, T], \quad \mathbf b = \hat P[T, M].$$

The matrix is block upper-triangular, so
$$\operatorname{eig}(\hat P) = \operatorname{eig}(Q) \cup \{1\}.$$

### 3.1 Why there is exactly one eigenvalue 1, and $|\lambda| < 1$ for the rest

It is enough to show that every eigenvalue of $Q$ has $|\lambda| < 1$. Two observations do this.

**$Q$ is block-triangular by groups.** Group the states by $(a, b)$. The only moves that leave a group lower $a$ or $b$, so with the groups in decreasing order $Q$ is block upper-triangular. Hence $\operatorname{eig}(Q)$ is the union of the eigenvalues of the diagonal blocks $B_{(a,b)}$.

**Each diagonal block is strictly substochastic in every row.** Inside a group only Re and C2 are possible, and neither is ever forbidden under either criterion. For a row of $B_{(a,b)}$,
$$\text{row sum} = \frac{k\rho + (a-k)(b-k)}{D}, \qquad D - \big(k\rho + (a-k)(b-k)\big) = \binom{a}{2} + \binom{b}{2} - \binom{k}{2} \ \ge\ \binom{b}{2} \ \ge\ 1.$$
The inequality uses $k \le a$ and that every transient allowed state has $a, b \ge 2$. So every row sum is strictly less than 1.

For a nonnegative matrix every eigenvalue satisfies $|\lambda| \le \max(\text{row sum})$. So every eigenvalue of every block, and hence of $Q$, has $|\lambda| < 1$.

**Consequences.**

- **1 is a simple eigenvalue.** $Q$ has no eigenvalue equal to 1, so the 1 comes only from $M$. This works because the forbidden states were removed. If they had been kept as absorbing states, each would add another eigenvalue 1.
- **No diagonalisation is needed.** Since $Q^m \to 0$,
$$\hat P^m = \begin{pmatrix} Q^m & (I + Q + \dots + Q^{m-1})\,\mathbf b \\ \mathbf 0^\top & 1 \end{pmatrix} \ \longrightarrow\ \begin{pmatrix} 0 & (I - Q)^{-1}\mathbf b \\ \mathbf 0^\top & 1 \end{pmatrix}.$$
  This holds whether or not $\hat P$ can be diagonalised. That matters, because the thesis's proof of Theorem 1 can diagonalise $P$ only when $R \ne 1$.
- **Probabilistic meaning.** $Q^m \to 0$ says the chain cannot stay among the allowed transient states forever: $a + b$ never increases, and from every state it has a positive chance of decreasing.

### 3.2 Why the eigenvalues come in pairs $\pm\beta$

**The parity argument.** Every move changes $a + b + k$ by an odd amount:

| Move | Re | C2 | left coalescence | right coalescence | C4 |
|---|---|---|---|---|---|
| change in $a+b+k$ | $-1$ | $+1$ | $-1$ | $-1$ | $-3$ |

So $Q$ only connects states of even parity with states of odd parity: the chain is bipartite.

Let $S = \operatorname{diag}\big((-1)^{a+b+k}\big)$. Every nonzero entry $Q_{xy}$ joins states of opposite parity, so $(SQS)_{xy} = -Q_{xy}$, i.e.
$$SQS = -Q.$$

Then $Q$ and $-Q$ are similar, so they have the same eigenvalues with the same multiplicities. Concretely, if $Q\mathbf v = \lambda\mathbf v$, then $Q(S\mathbf v) = -SQ\mathbf v = -\lambda\,(S\mathbf v)$, so $-\lambda$ is an eigenvalue with eigenvector $S\mathbf v$.

**Zero eigenvalues.** If the numbers of even and odd states differ, $Q$ has at least $|\#\text{even} - \#\text{odd}|$ zero eigenvalues. For $n = 2$ the three transient states $(2,2,2)$, $(2,2,1)$, $(2,2,0)$ are even, odd, even, which gives the spectrum $\{\beta, -\beta, 0\}$.

**Why 1 has no partner $-1$.** The absorbing self-loop at $M$ is not a move and does not change parity. It lies outside $Q$, so it breaks the symmetry only for the eigenvalue 1.

**The same thing at block level.** Each diagonal block $B_{(a,b)}$ contains only Re ($k-1$) and C2 ($k+1$) moves. It is tridiagonal with zero diagonal, so it is bipartite on the parity of $k$ and has its own $\pm$ pairs. Its off-diagonal products are positive, which makes the block similar to a symmetric matrix, so its eigenvalues are real.

For $n = 2$, $Q$ is the single block $B_{(2,2)}$:
$$Q = \begin{pmatrix} 0 & \frac{2R}{1+2R} & 0 \\ \frac{1}{3+R} & 0 & \frac{R}{3+R} \\ 0 & \frac23 & 0 \end{pmatrix}, \qquad \det(\lambda I - Q) = \lambda\left(\lambda^2 - \beta^2\right),$$
$$\beta^2 = \frac{2R}{(1+2R)(3+R)} + \frac{2R}{3(3+R)} = \frac{4R(R+2)}{3(3+R)(1+2R)},$$
which is the $\beta$ of Eq. (3.3).

**Example, $n = 4$, equal-heights criterion, $\rho = 0.8$.** Eigenvalues of each block:

| Group $(a,b)$ | Eigenvalues |
|---|---|
| (2,2) | $\pm 0.5498,\ 0$ |
| (2,3), (3,2) | $\pm 0.4158,\ 0$ |
| (2,4), (4,2) | $\pm 0.3130,\ 0$ |
| (3,3) | $\pm 0.4064,\ \pm 0.1254$ |
| (3,4), (4,3) | $\pm 0.3443,\ \pm 0.1096$ |
| (4,4) | $\pm 0.3337,\ \pm 0.1552,\ 0$ |

The largest, 0.5498, is the $n = 2$ value of $\beta$ at $\rho = 0.8$, coming from group (2,2).

---

## 4. The right eigenvector

### The left eigenvector

The last row of $\hat P$ is $\mathbf e_M^\top$, so $\mathbf e_M^\top \hat P = \mathbf e_M^\top$. In the long run, all mass that survives sits at the MRCA.

### Finding the right eigenvector

Solve $\hat P\,\mathbf r = \mathbf r$ with $\mathbf r = (\mathbf r_T, r_M)$:
$$Q\,\mathbf r_T + \mathbf b\, r_M = \mathbf r_T, \qquad r_M = r_M.$$
So $(I - Q)\,\mathbf r_T = \mathbf b\, r_M$. Since $I - Q$ is invertible (section 3.1),
$$\mathbf r_T = (I - Q)^{-1}\mathbf b\; r_M.$$
The eigenvector is unique up to scale, which reflects the fact that 1 is a simple eigenvalue. Set $r_M = 1$ and call the result $\mathbf h$.

### What $\mathbf h$ means

$h(x)$ is the probability of reaching the MRCA from state $x$ without entering a forbidden state. To see this, condition on the first step from $x$:

- a move to an allowed transient state $y$ (probability $Q(x,y)$) leads to success with probability $h(y)$;
- a move straight to the MRCA (probability $b(x)$) is a success;
- a move to a forbidden state is a failure, which contributes 0.

So
$$h(x) = \sum_{y \in T} Q(x,y)\,h(y) + b(x)\cdot 1,$$
which, row by row, is exactly $\hat P\mathbf h = \mathbf h$. In Markov-chain language, $\mathbf h$ is the harmonic function of absorption probabilities. **The right eigenvector for eigenvalue 1 is the vector of answers, one for each starting state.**

### Connection to the limit

With $\mathbf l = \mathbf e_M$ and $\mathbf l^\top\mathbf r = r_M = 1$,
$$\lim_{m\to\infty}\hat P^m = \frac{\mathbf r\,\mathbf l^\top}{\mathbf l^\top\mathbf r} = \mathbf h\,\mathbf e_M^\top.$$
This is the matrix whose last column is $\mathbf h$ and all other entries are 0, the same limit as in section 3.1. Row $x$ of the limit has total mass $h(x)$. The missing $1 - h(x)$ is the probability that the chain enters a forbidden state.

The required probability is $h(n,n,n)$, the (start, MRCA) entry of the limit.

### The fundamental matrix

Let $N = (I - Q)^{-1} = \sum_{m \ge 0} Q^m$. The entry $N(x,y)$ is the expected number of visits to $y$, starting from $x$, before the chain leaves $T$. Then
$$h(x) = \sum_y N(x,y)\,b(y),$$
i.e. a sum over the last allowed state before the MRCA.

Under both criteria the only allowed state that can step to the MRCA is $(2,2,2)$, by C4 with probability $1/(1+2\rho)$. So
$$P(\text{event}) = \frac{N\big((n,n,n),\,(2,2,2)\big)}{1 + 2\rho},$$
which is the expected number of visits to $(2,2,2)$ times the chance of finishing from there. I checked this numerically for $n = 2, 3, 4$.

### Worked example, $n = 2$

States 1, 2, 3 of the paper are $(2,2,2)$, $(2,2,1)$, $(2,2,0)$. The first-step equations are
$$h_1 = \tfrac{2R}{1+2R}h_2 + \tfrac{1}{1+2R}, \qquad h_2 = \tfrac{1}{3+R}h_1 + \tfrac{R}{3+R}h_3, \qquad h_3 = \tfrac23 h_2.$$

Substitute $h_3$ into the second equation: $h_2\left(3 + R - \tfrac{2R}{3}\right) = h_1$, so $h_2 = \dfrac{3h_1}{9+R}$.

Substitute this into the first equation: $h_1\big[(1+2R)(9+R) - 6R\big] = 9 + R$, and since $(1+2R)(9+R) - 6R = 2R^2 + 13R + 9$,
$$\mathbf h = \frac{(9+R,\ 3,\ 2,\ 2R^2 + 13R + 9)}{2R^2 + 13R + 9}, \qquad h_1 = \frac{R + 9}{2R^2 + 13R + 9}.$$
This is Theorem 2.

This vector is the first column of $\hat V$ in Eq. (3.9), up to scale. It also shows the typo there: the bottom-left entry of $\hat V$ should be $9 + 13R + 2R^2$, not $9 + 13R + R^2$.

### Computing it

Solving $(I - Q)\mathbf h = \mathbf b$ computes the right eigenvector directly. That is what `absorption_prob()`, `prob_equal_times()` and `prob_equal_heights()` in `transition_matrix.R` do. The linear solve is a sparse system with at most five nonzeros per row, so it needs no eigendecomposition and no matrix powers.
