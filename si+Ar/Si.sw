# Stillinger-Weber potential for Si and Ar (dummy)
# Original Si-Si-Si from F. H. Stillinger and T. A. Weber, PRB 31, 5262 (1985)
# Dummy parameters for ALL permutations of mixed triplets are set to zero.

# Format: element i, j, k, epsilon, sigma, a, lambda, gamma, cos(theta0), A, B, p, q, tol

# --- Real interaction for pure Si ---
Si  Si  Si  2.1686  2.0951  1.80  21.0  1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0

# --- Dummy entries for any triplet involving Ar (set Epsilon and Lambda to 0.0) ---

# Pure Ar
Ar  Ar  Ar  0.0     2.0951  1.80  0.0   1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0

# Two Ar, one Si (all permutations)
Ar  Ar  Si  0.0     2.0951  1.80  0.0   1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0
Ar  Si  Ar  0.0     2.0951  1.80  0.0   1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0
Si  Ar  Ar  0.0     2.0951  1.80  0.0   1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0

# One Ar, two Si (all permutations)
Ar  Si  Si  0.0     2.0951  1.80  0.0   1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0
Si  Ar  Si  0.0     2.0951  1.80  0.0   1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0
Si  Si  Ar  0.0     2.0951  1.80  0.0   1.20  -0.333333333333  7.049556277  0.6022245584  4.0  0.0  0.0