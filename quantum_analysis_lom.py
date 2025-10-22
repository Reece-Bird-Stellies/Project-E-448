import numpy as np
from scipy.constants import e, hbar, h
from typing import Tuple, Dict, Any, Optional
import pandas as pd

# =========================
# Ladder Operators
# =========================
def _a_op(N: int) -> np.ndarray:
    """
    Construct the bosonic annihilation operator for a truncated Fock space.
    
    The annihilation operator is defined as:
        a = sum_{n=1}^{N-1} sqrt(n) |n-1><n|
    
    Args:
        N: Dimension of the truncated Fock space (number of photon states)
    
    Returns:
        N×N complex matrix representing the annihilation operator
    """
    A = np.zeros((N, N), dtype=complex)
    for n in range(1, N):
        A[n-1, n] = np.sqrt(n)
    return A

def _adag_op(N: int) -> np.ndarray:
    """
    Construct the bosonic creation operator for a truncated Fock space.
    
    The creation operator is defined as:
        a† = sum_{n=0}^{N-2} sqrt(n+1) |n+1><n|
    
    Args:
        N: Dimension of the truncated Fock space (number of photon states)
    
    Returns:
        N×N complex matrix representing the creation operator
    """
    Adag = np.zeros((N, N), dtype=complex)
    for n in range(0, N-1):
        Adag[n+1, n] = np.sqrt(n+1)
    return Adag

# =========================
# System Parameters and Operators
# =========================
def _compute_constants(Lj, Lr, Cj, Cs, Cg, Cr, Nc, Nq):
    """
    Compute derived physical parameters and construct identity and ladder operators.
    
    This function calculates:
    - Total qubit capacitance (Csigma)
    - Capacitive coupling ratio (beta)
    - Resonator frequency (Wr)
    - RMS voltage fluctuations (Vrms)
    - Identity operators for cavity and qubit subspaces
    - Ladder operators for the cavity mode
    
    Args:
        Lj: Josephson inductance (H)
        Lr: Resonator inductance (H)
        Cj: Josephson junction capacitance (F)
        Cs: Shunt capacitance (F)
        Cg: Coupling capacitance (F)
        Cr: Resonator capacitance (F)
        Nc: Number of cavity Fock states to include
        Nq: Number of qubit energy levels to include
    
    Returns:
        Tuple containing:
        - Csigma: Total qubit capacitance (F)
        - beta: Capacitive coupling ratio (dimensionless)
        - Wr: Resonator angular frequency (rad/s)
        - Vrms: RMS voltage of resonator zero-point fluctuations (V)
        - Ic_hat: Identity operator for cavity subspace (Nc×Nc)
        - Iq_hat: Identity operator for qubit subspace (Nq×Nq)
        - a_hat: Cavity annihilation operator (Nc×Nc)
        - a_hat_dagger: Cavity creation operator (Nc×Nc)
    """
    Csigma = Cj + Cs + Cg
    beta = Cg / Csigma
    Wr = 1 / np.sqrt(Lr * Cr)
    Vrms = np.sqrt(hbar * Wr / (2 * Cr))
    
    Ic_hat = np.eye(Nc)
    Iq_hat = np.eye(Nq)
    a_hat = _a_op(Nc)
    a_hat_dagger = _adag_op(Nc)
    
    return (Csigma, beta, Wr, Vrms, Ic_hat, Iq_hat, a_hat, a_hat_dagger)

# =========================
# Cooper Pair Box Hamiltonian
# =========================
def _cpb_hamiltonian(Ec, Ej, ng, Nmax):
    """
    Construct the Cooper Pair Box (CPB) Hamiltonian in the charge basis.
    
    The CPB Hamiltonian is:
        H = 4*Ec*(n - ng)^2 - Ej*cos(φ)
    
    In the charge basis |N⟩ (N = -Nmax, ..., Nmax), this becomes:
        H = 4*Ec*sum_N (N - ng)^2 |N⟩⟨N| - (Ej/2)*sum_N (|N⟩⟨N+1| + |N+1⟩⟨N|)
    
    Args:
        Ec: Charging energy (J)
        Ej: Josephson energy (J)
        ng: Offset charge (dimensionless, in units of Cooper pairs)
        Nmax: Maximum charge number (basis spans -Nmax to +Nmax)
    
    Returns:
        Tuple containing:
        - H: CPB Hamiltonian matrix in charge basis (dim×dim, where dim = 2*Nmax+1)
        - n_charge: Charge number operator diagonal matrix
    """
    Nvals = np.arange(-Nmax, Nmax+1)
    dim = Nvals.size
    basis = np.eye(dim)
    
    n_charge = np.diag(Nvals.astype(float))
    
    # Charging energy term: 4*Ec*(n - ng)^2 (diagonal)
    proj_terms = []
    for i, N in enumerate(Nvals):
        weighted = (N - ng)**2 * np.outer(basis[i], basis[i])
        proj_terms.append(weighted)
    Hc = 4 * Ec * np.sum(proj_terms, axis=0)
    
    # Josephson energy term: -(Ej/2)*cos(φ) ≈ -(Ej/2)*(|N⟩⟨N+1| + h.c.) (off-diagonal)
    hop_terms = []
    for i in range(dim - 1):
        hop_terms.append(np.outer(basis[i], basis[i+1]) + np.outer(basis[i+1], basis[i]))
    Hj = -(Ej/2) * np.sum(hop_terms, axis=0)
    
    H = Hc + Hj
    return (H, n_charge)

# =========================
# Qubit Eigenstate Calculation
# =========================
def _extract_Wj_n_phi(Lj, Csigma, Nmax, Nq):
    """
    Diagonalize the CPB Hamiltonian and extract qubit energies and charge matrix elements.
    
    This function:
    1. Constructs the CPB Hamiltonian in the charge basis
    2. Diagonalizes to find energy eigenstates
    3. Keeps only the lowest Nq energy levels
    4. Transforms the charge operator to the energy eigenstate basis
    
    Args:
        Lj: Josephson inductance (H)
        Csigma: Total qubit capacitance (F)
        Nmax: Maximum charge number for CPB basis
        Nq: Number of lowest energy levels to keep
    
    Returns:
        Tuple containing:
        - Wj: Array of qubit energy eigenvalues for lowest Nq states (J)
        - n_phi: Charge operator matrix elements in energy eigenstate basis (Nq×Nq)
                 n_phi[i,j] = ⟨φ_i|n̂|φ_j⟩ where |φ_i⟩ are energy eigenstates
    """
    phi0 = hbar / (2*e)
    Ec = (e**2) / (2 * Csigma)
    Ej = phi0**2 / Lj
    ng = 0.0
    
    H, n_charge = _cpb_hamiltonian(Ec, Ej, ng, Nmax)
    E, U = np.linalg.eigh(H)
    idx = np.argsort(E)[:Nq]
    Wj = E[idx]
    Uq = U[:, idx]
    n_phi = Uq.conj().T @ n_charge @ Uq
    
    return (Wj, n_phi)

# =========================
# Coupling Strength Calculation
# =========================
def _compute_gij(Iq_hat, i, j, beta, Vrms, n_phi):
    """
    Compute the qubit-cavity coupling strength between energy levels i and j.
    
    The coupling strength in the dispersive/charge regime is:
        g_ij = (2*β*e*V_rms/ℏ) * ⟨φ_i|n̂|φ_j⟩
    
    where:
    - β is the capacitive coupling ratio Cg/(Cj+Cs+Cg)
    - V_rms is the RMS voltage of resonator zero-point fluctuations
    - ⟨φ_i|n̂|φ_j⟩ is the charge matrix element between eigenstates
    
    Args:
        Iq_hat: Identity operator for qubit subspace (unused, kept for consistency)
        i: Index of first energy level
        j: Index of second energy level
        beta: Capacitive coupling ratio (dimensionless)
        Vrms: RMS voltage fluctuations (V)
        n_phi: Charge operator matrix in eigenstate basis
    
    Returns:
        g_ij: Coupling strength between levels i and j (rad/s)
    """
    gij = (2*beta*e*Vrms / hbar) * n_phi[i, j]
    return gij

# =========================
# Total System Hamiltonian
# =========================
def compute_hamiltonian_lom(Lj, Lr, Cj, Cs, Cg, Cr, Nc, Nq, Nmax):
    """
    Construct the full Hamiltonian for a CPB qubit coupled to a resonator cavity.
    
    The total Hamiltonian in the dressed state basis is:
        H = H_r + H_q + H_i
    
    where:
    - H_r = ℏ*ω_r*a†a is the resonator Hamiltonian
    - H_q = sum_j E_j |j⟩⟨j| is the qubit Hamiltonian in energy eigenbasis
    - H_i = sum_ij ℏ*g_ij |i⟩⟨j| ⊗ (a + a†) is the interaction Hamiltonian
    
    The full Hilbert space has dimension Nq × Nc.
    Basis ordering: |qubit state⟩ ⊗ |cavity state⟩
    
    Args:
        Lj: Josephson inductance (H)
        Lr: Resonator inductance (H)
        Cj: Josephson junction capacitance (F)
        Cs: Shunt capacitance (F)
        Cg: Coupling capacitance (F)
        Cr: Resonator capacitance (F)
        Nc: Number of cavity Fock states
        Nq: Number of qubit energy levels
        Nmax: Maximum charge number for CPB diagonalization
    
    Returns:
        H: Total Hamiltonian matrix (Nq*Nc × Nq*Nc) in units of Joules
    """
    # Compute system parameters and operators
    (Csigma, beta, Wr, Vrms, Ic_hat, Iq_hat, a_hat, a_hat_dagger) = _compute_constants(
        Lj, Lr, Cj, Cs, Cg, Cr, Nc, Nq
    )
    
    # Extract qubit energies and charge matrix elements
    (Wj, n_phi) = _extract_Wj_n_phi(Lj, Csigma, Nmax, Nq)
    
    # Resonator Hamiltonian: H_r = ℏ*ω_r*a†a ⊗ I_q
    Hr = hbar * Wr * np.kron(Iq_hat, a_hat_dagger @ a_hat)
    
    # Qubit Hamiltonian: H_q = sum_j E_j |j⟩⟨j| ⊗ I_c
    proj_terms = []
    for j in range(Nq):
        x = Wj[j] * np.kron(np.outer(Iq_hat[j], Iq_hat[j]), Ic_hat)
        proj_terms.append(x)
    Hq = np.sum(proj_terms, axis=0)
    
    # Interaction Hamiltonian: H_i = sum_ij ℏ*g_ij |i⟩⟨j| ⊗ (a + a†)
    proj_terms2 = []
    for i in range(Nq):
        for j in range(Nq):
            gij = _compute_gij(Iq_hat, i, j, beta, Vrms, n_phi)
            x = gij * np.kron(np.outer(Iq_hat[i], Iq_hat[j]), a_hat + a_hat_dagger)
            proj_terms2.append(x)
    Hi = hbar * np.sum(proj_terms2, axis=0)
    
    # Total Hamiltonian
    H = Hr + Hq + Hi
    return H