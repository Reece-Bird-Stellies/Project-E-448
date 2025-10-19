from typing import Tuple, Dict, Any, Optional
import numpy as np
import pandas as pd
from math import factorial
from scipy.constants import e, hbar, h

PHI0_REDUCED: float = hbar / (2.0 * e)  # "reduced flux quantum" φ0 = ħ/2e (Wb)


"""
Two-mode (cavity + transmon) Hamiltonian builder, diagonalizer, and reporter.

- Uses participation ratios and Josephson inductance to form the junction phase operator.
- Linear part: ħ(ω_c n_c + ω_q n_q)
- Nonlinear part: -E_J [ cos(φ_J) + φ_J^2 / 2 ]
- Cosine computed via even-order Taylor series (stable for small φ_J as usual in transmon regimes).

Author: Reece Bird 
Credits: FILL IN
"""
# =========================
# Linear algebra utilities
# =========================
def ladder_ops(N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
def zero_point_phases(EJ_J: float, p_c: float, p_q: float, w_c: float, w_q: float) -> Tuple[float, float]:
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

def josephson_energy_from_Lj(Lj_H: float) -> float:
    """
    E_J = φ0^2 / L_J  with φ0 = ħ/(2e).
    Lj_H: Josephson inductance in Henry -> returns E_J in Joules.
    """
    if Lj_H <= 0:
        raise ValueError("Lj must be > 0.")
    
    return (PHI0_REDUCED ** 2) / Lj_H

def cosine_series(X: np.ndarray, order: int = 6) -> np.ndarray:
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
def build_two_mode_operators(N: int) -> Tuple[np.ndarray, ...]:
    """
    Build mode operators for a two-mode Hilbert space (cavity ⊗ qubit), each truncated to N.
    
    Returns a tuple in the order:
      (a_c, adag_c, n_c, a_q, adag_q, n_q)
      
    Shapes:
      - single-mode ops: (N,N)
      - two-mode ops: (N^2, N^2)
    """
    a, adag, n          = ladder_ops(N)
    I                   = np.eye(N, dtype=np.complex128)    # Same dimension truncation 
    Ic = I
    Iq = I
    a_c, adag_c, n_c    = np.kron(a, Iq), np.kron(adag, Iq), np.kron(n, Iq)
    a_q, adag_q, n_q    = np.kron(Ic, a), np.kron(Ic, adag), np.kron(Ic, n)

    return (a_c, adag_c, n_c, a_q, adag_q, n_q)


# =========================
# STEP 1: Hamiltonian builder
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
    EJ                                      = josephson_energy_from_Lj(Lj_H)

    a_c, adag_c, n_c, a_q, adag_q, n_q      = build_two_mode_operators(N)

    # Linear Hamiltonian (number operators)
    H_lin                                   = hbar * (w_c * n_c + w_q * n_q)

    # Nonlinear part: -EJ[cos(φJ) + φJ^2/2]
    phi_c, phi_q                            = zero_point_phases(EJ, p_c, p_q, w_c, w_q)
    phiJ                                    = phi_q * (a_q + adag_q) + phi_c * (a_c + adag_c)    #Change to phj_hat
    cos_phiJ                                = cosine_series(phiJ, order=cos_order)
    H_nl                                    = -EJ * (cos_phiJ + (phiJ @ phiJ)/2 )


    # Final Hamiltonian 
    H                                       = H_lin + H_nl

    return H


# =========================
# STEP 2: Diagonalization 
# =========================
def diagonalize(H: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Hermitian eigen-decomposition of H.
    Returns (eigenvalues ascending, eigenvectors columns).
    """
    E, V    = np.linalg.eigh(H) # eigh returns ascending already; keep explicit just in case
    idx     = np.argsort(E.real)
    E       = E[idx].real           # Can probably drop this
    V       = V[:, idx]
    return (E, V)


# =========================
# STEP 3: Helper 1
# =========================
def expectation(values: np.ndarray, op: np.ndarray) -> float:
    """
    ⟨ψ| op |ψ⟩ for a normalized state vector 'values'.
    """
    return float(np.vdot(values, op @ values).real)

# =========================
# STEP 3: Helper 2
# =========================
def characterize_eigenstates(
    V: np.ndarray,
    n_c: np.ndarray,
    n_q: np.ndarray,
    k: int = 6
) -> Dict[str, np.ndarray]:
    """
    Compute ⟨n_c⟩ and ⟨n_q⟩ for the lowest-k eigenstates.
    """
    k       = min(k, V.shape[1])
    nbar_c  = np.empty(k)
    nbar_q  = np.empty(k)
    for i in range(k):
        v           = V[:, i]
        nbar_c[i]   = expectation(v, n_c)
        nbar_q[i]   = expectation(v, n_q)
    return dict(nbar_c=nbar_c, nbar_q=nbar_q)


# =========================
# STEP 3: Summarize spectrum
# =========================
def summarize_spectrum(
    E_J: np.ndarray,
    E: np.ndarray,
    V: np.ndarray,
    n_c:float,
    n_q:float,
    k: int = 6
) -> Dict[str, Any]:
    """
    Extract core spectral features:
      - Qubit f01 (GHz) and anharmonicity α_q (MHz) from first spacings.
      - Heuristic cavity 1-photon line and its anharmonicity α_c (MHz) via ⟨n_c⟩.
    """

    k               = min(k, len(E))
    dE              = np.diff(E[:max(3, k)])  # need at least E0,E1,E2 for α_q   # Can probably remove this

    # Qubit f01, α_q
    w01             = (E[1] - E[0]) / h           # rad/s
    w12             = (E[2] - E[1]) / h
    f01_GHz         = w01 /  1e9
    alpha_q_MHz     = (w12 - w01)  / 1e6

    # Mode characters via ⟨n_c⟩
    chars           = characterize_eigenstates(V, n_c, n_q, k=min(20, len(E)))
    nbar_c          = chars["nbar_c"]

    # Pick a "cavity-like 1-photon" as state with ⟨n_c⟩ in ~[0.4,1.6]
    i1              = next((i for i, val in enumerate(nbar_c) if 0.4 <= val <= 1.6), None)
    # Pick a "cavity-like 2-photon" as ⟨n_c⟩ in ~[1.6,2.6]
    i2              = next((i for i, val in enumerate(nbar_c) if i != i1 and 1.6 <= val <= 2.6), None)

    f_cav_spec_GHz  = None
    alpha_c_MHz     = None
    if i1 is not None:
        f_cav_spec_GHz = (E[i1] - E[0]) / h / 1e9
        if i2 is not None:
            f12     = (E[i2] - E[i1]) / h / 1e6  # MHz
            f01c    = (E[i1] - E[0]) / h / 1e6  # MHz
            alpha_c_MHz = f12 - f01c

    return dict(
        f01_GHz         = f01_GHz,
        alpha_q_MHz     = alpha_q_MHz,
        f_cav_spec_GHz  = f_cav_spec_GHz,
        alpha_c_MHz     = (alpha_c_MHz if alpha_c_MHz is not None else 0.0),
        dE_first        = dE                    # Can probably remove this 
    )


# =========================
# STEP 4: Compute derived scalars
# =========================
def compute_scalars(
    Hdict: Dict[str, Any],
    spectrum: Dict[str, Any],
    Lj_H: Optional[float] = None,
    g_RWA_MHz: Optional[float] = None
) -> Dict[str, float]:
    """
    Compute derived quantities used in reporting.
    """
    EJ              = float(Hdict["EJ"])
    w_c, w_q        = float(Hdict["w_c"]), float(Hdict["w_q"])
    f_c_lin_GHz     = w_c / (2.0 * np.pi) / 1e9
    f_q_lin_GHz     = w_q / (2.0 * np.pi) / 1e9
    phi_c, phi_q    = float(Hdict["phi_c"]), float(Hdict["phi_q"])

    # Participations recovered from ZPF (sanity check / echo)
    p_c             = (phi_c ** 2) * (2.0 * EJ) / (hbar * w_c)
    p_q             = (phi_q ** 2) * (2.0 * EJ) / (hbar * w_q)
    p_sum           = p_c + p_q
    p_c_norm        = p_c / p_sum if p_sum > 0 else np.nan
    p_q_norm        = p_q / p_sum if p_sum > 0 else np.nan

    # Spectral results
    f_q_spec_GHz    = float(spectrum["f01_GHz"])
    alpha_q_MHz     = float(spectrum["alpha_q_MHz"])
    f_c_spec_GHz    = float(spectrum["f_cav_spec_GHz"]) if spectrum["f_cav_spec_GHz"] is not None else f_c_lin_GHz
    alpha_c_MHz     = float(spectrum["alpha_c_MHz"])

    # Detuning (spectral) in GHz
    detuning_GHz    = f_c_lin_GHz - f_q_spec_GHz

    # Flux ZPF at junction (ground)
    flux_zpf_sq     = phi_c ** 2 + phi_q ** 2

    # Transmon estimates: Ec/h ≈ -α_q
    Ec_GHz          = -alpha_q_MHz / 1e3
    EJ_over_h_GHz   = EJ / h / 1e9
    EJ_over_EC      = EJ_over_h_GHz / Ec_GHz if Ec_GHz > 0 else np.nan

    # C_qubit from Ec:  Ec = e^2/(2 Cq) ⇒ Cq = e^2/(2 h Ec)  (Ec in Hz)
    Ec_Hz           = Ec_GHz * 1e9
    Cq_F            = (e ** 2) / (2.0 * h * Ec_Hz) if Ec_Hz > 0 else np.nan
    Cq_fF           = Cq_F * 1e15 if np.isfinite(Cq_F) else np.nan

    return dict(
        Ej_GHz                          = EJ_over_h_GHz,
        Ec_GHz                          = Ec_GHz,
        Ej_over_Ec                      = EJ_over_EC,
        Lj_nH                           = (Lj_H * 1e9 if Lj_H is not None else np.nan),
        C_qubit_fF                      = Cq_fF,   
        Participation_Qubit             = p_q,
        Participation_Qubit_Normalized  = p_q_norm,
        Participation_Cavity            = p_c,
        Participation_Cavity_Normalized = p_c_norm,
        Linear_Qubit_Frequency_GHz      = f_q_lin_GHz,
        Qubit_Frequency_GHz             = f_q_spec_GHz,
        Qubit_Anharmonicity_MHz         = alpha_q_MHz,
        Linear_Cavity_Frequency_GHz     = f_c_lin_GHz,
        Cavity_Frequency_GHz            = f_c_spec_GHz,
        Cavity_Anharmonicity_MHz        = alpha_c_MHz,
        Coupling_g_RWA_MHz              = (float(g_RWA_MHz) if g_RWA_MHz is not None else np.nan),
        Detuning_GHz                    = detuning_GHz,
        Flux_ZPF_squared                = flux_zpf_sq
    )


# =========================
# STEP 5: Return Hdict, spectrum, and a report DataFrame
# =========================
def make_report_dataframe(scalars: Dict[str, float]) -> pd.DataFrame:
    """Tabulate the key outputs."""
    rows = [
        ("Ej",                               scalars["Ej_GHz"],                         "GHz"),
        ("Ec",                               scalars["Ec_GHz"],                         "GHz"),
        ("Ej/Ec",                            scalars["Ej_over_Ec"],                     "unitless"),
        ("Lj",                               scalars["Lj_nH"],                          "nH"),
        ("C_qubit",                          scalars["C_qubit_fF"],                     "fF"),
        ("Participation Ratio Qubit",        scalars["Participation_Qubit"],            "unitless"),
        ("Participation Ratio Qubit (Normalized)", scalars["Participation_Qubit_Normalized"], "unitless"),
        ("Participation Ratio Cavity",       scalars["Participation_Cavity"],           "unitless"),
        ("Participation Ratio Cavity (Normalized)", scalars["Participation_Cavity_Normalized"], "unitless"),
        ("Linear Qubit Frequency",           scalars["Linear_Qubit_Frequency_GHz"],     "GHz"),
        ("Qubit Frequency",                  scalars["Qubit_Frequency_GHz"],            "GHz"),
        ("Qubit Anharmonicity",              scalars["Qubit_Anharmonicity_MHz"],        "MHz"),
        ("Linear Cavity Frequency",          scalars["Linear_Cavity_Frequency_GHz"],    "GHz"),
        ("Cavity Frequency",                 scalars["Cavity_Frequency_GHz"],           "GHz"),
        ("Cavity Anharmonicity",             scalars["Cavity_Anharmonicity_MHz"],       "MHz"),
        ("Cavity-Qubit Coupling [RWA], g",   scalars["Coupling_g_RWA_MHz"],             "MHz"),
        ("Detuning",                         scalars["Detuning_GHz"],                   "GHz"),
        ("Flux_ZPF_squared",                 scalars["Flux_ZPF_squared"],               "unitless"),
    ]
    return pd.DataFrame(rows, columns=["Parameter", "Value", "Unit"])


# =========================
# Orchestration helper
# =========================
def solve_and_report(
    Lj_H: float,
    pc: float,
    pq: float,
    fc_Hz: float,
    fq_Hz: float,
    N: int = 10,
    cos_order: int = 6,
    k_levels: int = 6,
    g_RWA_MHz: Optional[float] = None
) -> Tuple[Dict[str, Any], Dict[str, Any], pd.DataFrame]:
    """
    End-to-end:
      1) Build Hamiltonian
      2) Diagonalize
      3) Summarize spectrum
      4) Compute derived scalars
      5) Return Hdict, spectrum, and a report DataFrame
    """
    Hdict       = compute_hamiltonian_epr(Lj_H, pc, pq, fc_Hz, fq_Hz, N=N, cos_order=cos_order)
    E, V        = diagonalize(Hdict["H"])
    spectrum    = summarize_spectrum(Hdict["EJ"], E, V, n_c=Hdict['n_c'], n_q=Hdict["n_q"], k=k_levels)
    scalars     = compute_scalars(Hdict, spectrum, Lj_H=Lj_H, g_RWA_MHz=g_RWA_MHz)
    df          = make_report_dataframe(scalars)
    return (Hdict, spectrum, df)