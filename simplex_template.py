import numpy as np
import sys
import argparse

eps = 0.00001

def lu_factor(A):
    A = A.copy().astype(float)
    n = A.shape[0]
    piv = np.arange(n)
    for k in range(n - 1):
        p = k + np.argmax(np.abs(A[k:, k]))
        if abs(A[p, k]) < eps:
            continue
        if p != k:
            A[[k, p]] = A[[p, k]]
            piv[[k, p]] = piv[[p, k]]
        A[k+1:, k] /= A[k, k]
        A[k+1:, k+1:] -= np.outer(A[k+1:, k], A[k, k+1:])
    return A, piv

def lu_solve(LU, piv, b):
    b = b.copy().astype(float)
    b = b[piv]
    n = LU.shape[0]

    for i in range(n):
        b[i+1:] -= LU[i+1:, i] * b[i]

    for i in range(n - 1, -1, -1):
        if abs(LU[i, i]) < eps:
            b[i] = b[i] / (LU[i, i] if abs(LU[i, i]) > eps else eps)
        else:
            b[i] /= LU[i, i]
        b[:i] -= LU[:i, i] * b[i]
    return b

def basis_factor(A, basis):
    B = A[:, basis]
    return lu_factor(B)

def basis_replace_column(A, basis, LU, piv, leave_pos, entering):
    basis[leave_pos] = entering
    LU, piv = basis_factor(A, basis)
    return basis, LU, piv

def PrimalSimplex(c, A, b, basis=None, nbasis=None):
    m, n = A.shape
    n -= m
    if basis is None or nbasis is None:
        basis = list(range(n, n + m))
        nbasis = list(range(0, n))

    LU, piv = basis_factor(A, basis)

    while True:
        cB = c[basis]
        B = A[:, basis]

        LUt, pivt = lu_factor(B.T)
        pi = lu_solve(LUt, pivt, cB)

        N = A[:, nbasis]
        reduced_cost = c[nbasis] - pi @ N

        entering_index = None
        for idx, col in enumerate(nbasis):
            if reduced_cost[idx] > eps:
                entering_index = col
                enter_pos = idx
                break

        if entering_index is None:
            x = np.zeros(n + m)
            xB = lu_solve(LU, piv, b)
            x[basis] = xB
            return "optimal", x[:n], float(c @ x)

        aj = A[:, entering_index]
        d = lu_solve(LU, piv, aj)

        if np.all(d <= eps):
            return "unbounded", None, None

        xB = lu_solve(LU, piv, b)
        ratios = np.full_like(d, np.inf, dtype=float)
        np.divide(xB, d, out=ratios, where=d > eps)

        min_ratio = np.min(ratios)
        cands = np.where(np.abs(ratios - min_ratio) <= eps)[0]
        chosen_pos = min(cands, key=lambda p: basis[p])

        leaving_index = basis[chosen_pos]

        basis, LU, piv = basis_replace_column(A, basis, LU, piv, chosen_pos, entering_index)
        nbasis[enter_pos] = leaving_index

def Phase1(c, A, b):
    m, n = A.shape

    new_c = np.concatenate([np.zeros(n), np.zeros(m), -np.ones(1)])
    new_A = np.concatenate([A, np.eye(m), -np.ones((m, 1))], axis=1)

    nm_sum = n + m

    basis  = list(range(n, n + m))
    nbasis = [*range(n), nm_sum]

    j = int(np.argmin(b))
    basis[j] = nm_sum
    nbasis[-1] = n + j

    status, x_aux, obj = PrimalSimplex(new_c, new_A, b, basis.copy(), nbasis.copy())
    if status != "optimal" or obj < -eps:
        return "infeasible", None, None

    if nm_sum in basis:
        r = basis.index(nm_sum)
        for j in nbasis:
            if abs(new_A[r, j]) > eps:
                basis[r] = j
                nbasis.remove(j)
                break
        else:
            nbasis.remove(nm_sum)
    else:
        nbasis = [j for j in nbasis if j != nm_sum]

    c = np.concatenate([c, np.zeros(m)])
    A = np.concatenate([A, np.eye(m)], axis=1)

    LU, piv = basis_factor(A, basis)
    xB = lu_solve(LU, piv, b)
    if np.any(xB < -eps):
        return "infeasible", None, None

    return PrimalSimplex(c, A, b, basis, nbasis)

def DualSimplex(c, A, b, basis=None, nbasis=None):
    m, n = A.shape
    n -= m

    if basis is None or nbasis is None:
        basis = list(range(n, n + m))
        nbasis = list(range(n))

    B_inv = np.linalg.inv(A[:, basis])
    N = A[:, nbasis]

    almost_inf = 10**10

    while True:
        w = B_inv.T @ c[basis]
        rN = c[nbasis] - (N.T @ w)

        xN = np.where(rN > eps, almost_inf, 0.0)
        d = np.where(xN == almost_inf, -1.0, 1.0)

        xB = B_inv @ (b - N @ xN)
        r = d * rN
        if np.all(xB >= -eps):
            break

        candidate = np.argmin(xB)
        a = B_inv @ N
        d_r = - a[candidate, :] * d
        valid = d_r > eps
        if not np.any(valid):
            return "infeasible", None, None

        theta = -r[valid] / d_r[valid]
        j = np.flatnonzero(valid)[np.argmin(theta)]

        basis[candidate], nbasis[j] = nbasis[j], basis[candidate]

        B = A[:, basis]
        N = A[:, nbasis]

        u = B[:, candidate].reshape(-1, 1) - N[:, j].reshape(-1, 1)
        v = np.zeros((m, 1)); v[candidate, 0] = 1
        B_inv -= (B_inv @ u) @ (v.T @ B_inv) / (1 + (v.T @ B_inv @ u)[0, 0])

    x = np.zeros(n+m)
    x[basis] = xB
    x[nbasis] = xN

    sol = x[:n]
    z = float((c[basis] @ xB))

    if np.any(sol == almost_inf):
        return "unbounded", None, None

    return "optimal", sol, z


def Solve(c, A, b, mode="dual"):
    m = A.shape[0]
    c = np.concatenate([c, np.zeros(m)])
    A = np.concatenate([A, np.eye(m)], axis=1)
    if mode == "primal" or (mode == "auto" and np.all(b >= -eps)):
        if np.all(b >= -eps):
            return PrimalSimplex(c, A, b)
        else:
            return Phase1(c, A, b)
    elif mode == "dual" or (mode == "auto" and np.all(b >= -eps)):
        return DualSimplex(c, A, b)
    else:
        raise ValueError("Invalid mode: mode must be one of 'primal', 'dual' or 'auto'")

def proc_cmd():
    parser = argparse.ArgumentParser(description="Solve a linear program using the Primal Simplex method.")
    parser.add_argument("filename", type=str, help="Input file containing the LP problem.")
    return parser.parse_args()

def main():
    args = proc_cmd()
    with open(args.filename, 'r', encoding='utf-8') as f:
        n, m = map(int, f.readline().split())
        c = np.array(list(map(float, f.readline().split())))
        A = []
        b = []
        for _ in range(m):
            *row, bi = map(float, f.readline().split())
            A.append(row)
            b.append(bi)
        A = np.array(A)
        b = np.array(b)
    print("n =", n)
    print("m =", m)
    print("c =", c)
    print("A =", A)
    print("b =", b)

    print("Solving the linear program using the Primal Simplex method...\n")
    status, solution, objective = Solve(c, A, b)
    print("\nResult:")
    print("Status:", status)
    if status == "optimal":
        print("Optimal solution x* =", solution)
        print("Optimal value =", objective)
    elif status == "unbounded":
        print("The problem is unbounded.")
    else:
        print("No solution found.")

if __name__ == '__main__':
    main()
