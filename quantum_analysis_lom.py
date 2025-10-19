import numpy as np
from scipy.constants import e, hbar, h
from typing import Tuple, Dict, Any, Optional
import pandas as pd

# =========================
# Constants helper
# =========================
def _a_op(N: int) -> np.ndarray:
    """
    Construct the annihilation operator a_N using explicit loops.
    a_N = sum_{n=1}^{N-1} sqrt(n) |n-1><n|
    """
    A = np.zeros((N, N), dtype=complex)
    for n in range(1, N):  # 1 → N-1
        A[n-1, n] = np.sqrt(n)
    return A

def _adag_op(N: int) -> np.ndarray:
    """
    Construct the creation operator a_N^† using explicit loops.
    a_N^† = sum_{n=0}^{N-2} sqrt(n+1) |n+1><n|
    """
    Adag = np.zeros((N, N), dtype=complex)
    for n in range(0, N-1):  # 0 → N-2
        Adag[n+1, n] = np.sqrt(n+1)
    return Adag

# =========================
# Constants
# =========================
def _compute_constants(Lj,Lr,Cj,Cs,Cg,Cr,Nc,Nq):
    Csigma          = Cj + Cs + Cg
    beta            = Cg / Csigma
    Wr              = 1 / np.sqrt(Lr * Cr)
    Vrms            = np.sqrt(hbar * Wr / (2 * Cr))
    #a_hat2          = destroy(Nc)
    #a_hat_dagger2   = a_hat.dag()
    Ic_hat          = np.eye(Nc)
    Iq_hat          = np.eye(Nq)
    a_hat           = _a_op(Nc)
    a_hat_dagger    = _adag_op(Nc)
    return (Csigma,beta, Wr, Vrms, Ic_hat, Iq_hat, a_hat, a_hat_dagger)

# =========================
# extract_Wj_n_hat helper
# =========================
def _cpb_hamiltonian(Ec, Ej, ng, Nmax):
    """
    Build H = 4 Ec (n - ng)^2 - Ej cos(phi) in the charge basis |N>, N=-Nmax..Nmax,
    using outer products and np.sum for the summations.
    Returns: H (dim x dim matrix), Nvals (basis labels)
    """
    Nvals           = np.arange(-Nmax, Nmax+1)           # N = -Nmax,...,Nmax
    dim             = Nvals.size
    basis           = np.eye(dim)                           # |N_i> are canonical unit vectors

    n_charge        = np.diag(Nvals.astype(float))

    # ---- Charging term: 4 Ec sum_N (N - ng)^2 |N><N|  (diagonal)
    proj_terms      = []       # list to collect all projectors
    for i, N in enumerate(Nvals):
        weighted    =  (N - ng)**2 * np.outer(basis[i], basis[i])# weighted term: (N - ng)^2 |N><N|
        proj_terms.append(weighted)
    Hc              = 4 * Ec * np.sum(proj_terms, axis=0)

    # ---- Josephson term: -(Ej/2) sum_N ( |N><N+1| + |N+1><N| ) (off-diagonal)
    hop_terms       = []
    for i in range(dim - 1):
        hop_terms.append(np.outer(basis[i],   basis[i+1]) +  np.outer(basis[i+1], basis[i]))  # |N><N+1| + |N+1><N|
    Hj              = -(Ej/2) * np.sum(hop_terms, axis=0)
    H = Hc + Hj
    return (H, n_charge)

# =========================
# extract_Wj_n_hat 
# =========================
def _extract_Wj_n_phi(Lj, Csigma, Nmax, Nq):
    phi0 = hbar/(2*e)
    Ec   = (e**2) / (2 * Csigma)
    Ej   = phi0**2 / Lj
    ng   = 0.0

    H, n_charge = _cpb_hamiltonian(Ec, Ej, ng, Nmax)
    E, U = np.linalg.eigh(H)                # columns U[:,j] = |φ_j>
    idx  = np.argsort(E)[:Nq]               # keep lowest Nq levels
    Wj   = E[idx]                           # energies (J)
    Uq   = U[:, idx]                        # eigenvectors for those levels
    n_phi = Uq.conj().T @ n_charge @ Uq     # <φ_i| n̂ |φ_j> (Nq×Nq)   #NOT CONVINCED THIS RIGHT
    return (Wj, n_phi)

# =========================
# compute_gij 
# =========================
def _compute_gij(Iq_hat, i, j, beta, Vrms, n_phi):
    gij = (2*beta*e*Vrms / hbar) * n_phi[i, j]
    return gij

# =========================
# Compute Hamiltonian 
# =========================
def compute_hamiltonian_lom(Lj, Lr, Cj, Cs, Cg, Cr, Nc, Nq, Nmax):
    (Csigma, beta, Wr, Vrms, Ic_hat, Iq_hat, a_hat, a_hat_dagger) = _compute_constants(
        Lj, Lr, Cj, Cs, Cg, Cr, Nc, Nq
    )
    (Wj, n_phi) = _extract_Wj_n_phi(Lj, Csigma, Nmax, Nq)  # energies (J) and n_phi (Nq×Nq)

    # Resonator Hamiltonian
    Hr = hbar * Wr * np.kron(Iq_hat, a_hat_dagger @ a_hat)

    # Qubit Hamiltonian
    proj_terms = []
    for j in range(Nq):
        x = Wj[j] * np.kron(np.outer(Iq_hat[j], Iq_hat[j]), Ic_hat)  # energies => no extra ħ
        proj_terms.append(x)
    Hq = np.sum(proj_terms, axis=0)

    # Interaction Hamiltonian
    proj_terms2 = []
    for i in range(Nq):
        for j in range(Nq):
            gij = _compute_gij(Iq_hat, i, j, beta, Vrms, n_phi)  # uses n_phi directly
            x = gij * np.kron(np.outer(Iq_hat[i], Iq_hat[j]), a_hat + a_hat_dagger)
            proj_terms2.append(x)
    Hi = hbar * np.sum(proj_terms2, axis=0)  # multiply by ħ to get energy units

    # Total Hamiltonian
    H = Hr + Hq + Hi
    return H


