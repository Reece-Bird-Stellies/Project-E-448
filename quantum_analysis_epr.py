from typing import Tuple, Dict, Any, Optional
import numpy as np
from math import factorial
from scipy.constants import e, hbar, h

PHI0_REDUCED: float = hbar / (2.0 * e)  # "reduced flux quantum" φ0 = ħ/2e (Wb)

# =========================
# Linear algebra utilities
# =========================
def _ladder_ops(N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Create annihilation (a), creation (adag), and number (n = a†a) operators
    in the |n> Fock basis of dimension N.
    """
    if N < 2:
        raise ValueError("Truncation N must be >= 2.")
    
    a               = np.zeros((N, N), dtype=np.complex128)
    for n in range(1, N):   
        a[n - 1, n] = np.sqrt(n)  # <n-1| a |n> = sqrt(n)
    adag            = a.conj().T
    n_op            = adag @ a
    return (a, adag, n_op)                # Includes classsic floating point error. idk if this is problem (15.000000000000002+0j) != 15


# =========================
# Physics helpers
# =========================
def _zero_point_phases(EJ_J: float, p_c: float, p_q: float, w_c: float, w_q: float) -> Tuple[float, float]:
    """
    φ_c^2 = p_c ħ ω_c / (2 E_J), φ_q^2 = p_q ħ ω_q / (2 E_J).
    Returns (φ_c, φ_q), dimensionless phase amplitudes at the junction.
    """
    if EJ_J <= 0:
        raise ValueError("EJ must be > 0.")
    if any(x < 0 for x in (p_c, p_q)):
        raise ValueError("Participations must be >= 0.")
    if any(x <= 0 for x in (w_c, w_q)):
        raise ValueError("Frequencies (rad/s) must be > 0.")

    phi_c = np.sqrt(p_c * hbar * w_c / (2.0 * EJ_J))
    phi_q = np.sqrt(p_q * hbar * w_q / (2.0 * EJ_J))
    return (float(phi_c), float(phi_q))

def _josephson_energy_from_lj(Lj_H: float) -> float:
    """
    E_J = φ0^2 / L_J  with φ0 = ħ/(2e).
    Lj_H: Josephson inductance in Henry -> returns E_J in Joules.
    """
    if Lj_H <= 0:
        raise ValueError("Lj must be > 0.")
    
    return (PHI0_REDUCED ** 2) / Lj_H

def _cosine_series(X: np.ndarray, order: int = 6) -> np.ndarray:
    """
    cos(X) via even-order Taylor series up to 'order' (which will be clamped to the nearest even ≤ order).
    cos(X) = I - X^2/2! + X^4/4! - X^6/6! + ...

    Notes:
      - For transmon regimes (small φ_J ZPF) low even orders are typically adequate.
      - If you need stronger guarantees for larger norms, consider scaling & squaring with Pade; kept explicit by request.
    """
    if order < 0:
        raise ValueError("order must be non-negative.")
    if X.ndim != 2 or X.shape[0] != X.shape[1]:
        raise ValueError("X must be a square matrix.")

    # Force even order: 0, 2, 4, 6, ...
    kmax    = (order // 2) * 2
    I       = np.eye(X.shape[0], dtype=X.dtype)
    
    if kmax == 0:
        return I

    X2      = X @ X
    cosX    = I.copy()
    power   = X2.copy()   # X^2
    sign    = -1
    steps   = kmax // 2
    for m in range(1, steps + 1):
        cosX = cosX + (sign / factorial(2 * m)) * power
        if m < steps:
            power = power @ X2  # X^{2(m+1)}
        sign *= -1
    return cosX


# =========================
# Operator construction
# =========================
def _build_two_mode_operators(N: int) -> Tuple[np.ndarray, ...]:
    """
    Build mode operators for a two-mode Hilbert space (cavity ⊗ qubit), each truncated to N.
    
    Returns a tuple in the order:
      (a_c, adag_c, n_c, a_q, adag_q, n_q)
      
    Shapes:
      - single-mode ops: (N,N)
      - two-mode ops: (N^2, N^2)
    """
    a, adag, n          = _ladder_ops(N)
    I                   = np.eye(N, dtype=np.complex128)    # Same dimension truncation 
    Ic = I
    Iq = I
    a_c, adag_c, n_c    = np.kron(a, Iq), np.kron(adag, Iq), np.kron(n, Iq)
    a_q, adag_q, n_q    = np.kron(Ic, a), np.kron(Ic, adag), np.kron(Ic, n)

    return (a_c, adag_c, n_c, a_q, adag_q, n_q)


# =========================
# Hamiltonian builder
# =========================
def compute_hamiltonian_epr(
    Lj_H: float,
    p_c: float,
    p_q: float,
    f_c_Hz: float,
    f_q_Hz: float,
    N: int = 10,
    cos_order: int = 6
) -> Dict[str, Any]:
    """
    Build H = H_lin + H_nl for the coupled cavity-transmon system.

    Inputs:
      Lj_H       : Josephson inductance [H]
      p_c, p_q   : participation ratios (dimensionless)
      f_c_Hz, f_q_Hz : *linear* frequencies [Hz]
      N          : truncation per mode
      cos_order  : even max order for Taylor cosine

    Returns dict with:
      H, H_lin, H_nl [J]
      operators (a_c, adag_c, n_c, a_q, adag_q, n_q)
      scalars: EJ [J], w_c, w_q [rad/s], phi_c, phi_q (dimensionless), N, cos_order
    """
    if any(x <= 0 for x in (f_c_Hz, f_q_Hz)):
        raise ValueError("Frequencies (Hz) must be > 0.")

    w_c                                     = 2.0 * np.pi * f_c_Hz
    w_q                                     = 2.0 * np.pi * f_q_Hz
    EJ                                      = _josephson_energy_from_lj(Lj_H)

    a_c, adag_c, n_c, a_q, adag_q, n_q      = _build_two_mode_operators(N)

    # Linear Hamiltonian (number operators)
    H_lin                                   = hbar * (w_c * n_c + w_q * n_q)

    # Nonlinear part: -EJ[cos(φJ) + φJ^2/2]
    phi_c, phi_q                            = _zero_point_phases(EJ, p_c, p_q, w_c, w_q)
    phiJ                                    = phi_q * (a_q + adag_q) + phi_c * (a_c + adag_c)    #Change to phj_hat
    cos_phiJ                                = _cosine_series(phiJ, order=cos_order)
    H_nl                                    = -EJ * (cos_phiJ + (phiJ @ phiJ)/2 )


    # Final Hamiltonian 
    H                                       = H_lin + H_nl

    return H

