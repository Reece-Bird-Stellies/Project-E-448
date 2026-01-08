import numpy as np
from scipy.constants import e, hbar

# ============================================================
# Ladder Operators (Truncated Bosonic Fock Space)
# ============================================================
def _a_op(N: int) -> np.ndarray:
    """
    Construct the bosonic annihilation operator â for a truncated Fock space.

    The annihilation operator is defined as:
            N-1
        â =  Σ √n |n-1⟩⟨n|
            n=1

    Args:
        N : Dimension of the truncated Fock space

    Returns:
        A : NxN complex matrix representing â
    """
    A = np.zeros((N, N), dtype=complex) # Allocate NxN zero matrix
    for n in range(1, N):
        A[n-1, n] = np.sqrt(n)          # Populate superdiagonal with √n
    return A


def _adag_op(N: int) -> np.ndarray:
    """
    Construct the bosonic creation operator â† for a truncated Fock space.

    The creation operator is defined as:
            N-2
        â† = Σ √(n+1) |n+1⟩⟨n|
            n=0
    Args:
        N : Dimension of the truncated Fock space

    Returns:
        A_dag : NxN complex matrix representing â†
    """
    A_dag = np.zeros((N, N), dtype=complex)
    for n in range(0, N-1):
        A_dag[n+1, n] = np.sqrt(n+1)    # Populate subdiagonal with √(n+1)
    return A_dag


# ============================================================
# System Constants and Basic Operators
# ============================================================
def _compute_constants(Lr, Cj, Cs, Cg, Cr, Nc, Nq):
    """
    Compute derived circuit parameters and construct basic operators.

    This routine computes:
    - Total qubit capacitance ΣC
    - Capacitive coupling ratio β
    - Resonator angular frequency ω_r
    - Zero-point voltage fluctuations V_zpf
    - Identity operators for cavity and qubit subspaces
    - Cavity ladder operators â and â†

    Args:
        Lr : Resonator inductance (H)
        Cj : Junction capacitance (F)
        Cs : Shunt capacitance (F)
        Cg : Coupling capacitance (F)
        Cr : Resonator capacitance (F)
        Nc : Number of cavity Fock states
        Nq : Number of qubit energy levels

    Returns:
        Csigma       : Total qubit capacitance ΣC
        beta         : Capacitive coupling ratio β
        Wr           : Resonator angular frequency ω_r (rad/s)
        Vzpf         : RMS zero-point voltage fluctuations V_zpf
        Ic_hat       : Cavity identity operator
        Iq_hat       : Qubit identity operator
        a_hat        : Cavity annihilation operator â
        a_hat_dagger : Cavity creation operator â†
    """
    Csigma       = Cj + Cs + Cg
    beta         = Cg / Csigma
    Wr           = 1 / np.sqrt(Lr * Cr)
    Vzpf         = np.sqrt((hbar * Wr) / (2 * Cr))

    Ic_hat       = np.eye(Nc)  # NcxNc identity matrix
    Iq_hat       = np.eye(Nq)  # NqxNq identity matrix
    a_hat        = _a_op(Nc)
    a_hat_dagger = _adag_op(Nc)

    return (Csigma, beta, Wr, Vzpf, Ic_hat, Iq_hat, a_hat, a_hat_dagger)


# ============================================================
# Cooper Pair Box Hamiltonian (Charge Basis)
# ============================================================
def _cpb_hamiltonian(Ec, Ej, ng, Nmax):
    """
    Construct the Cooper Pair Box (CPB) Hamiltonian in the charge basis.

    The CPB Hamiltonian is:
        Ĥ = 4E_c (n̂ - n_g)^2 - E_J cos(φ̂ )

    In the truncated charge basis |N⟩ with N ∈ [-Nmax, …, Nmax],
    this becomes:

        Charging term:
            4E_c (N - n_g)^2 |N⟩⟨N|

        Josephson term:
            -(E_J / 2) ( |N⟩⟨N+1| + |N+1⟩⟨N| )

    Args:
        Ec   : Charging energy E_c (J)
        Ej   : Josephson energy E_J (J)
        ng   : Offset charge n_g (dimensionless)
        Nmax : Maximum charge index

    Returns:
        H        : CPB Hamiltonian matrix
        n_charge : Charge number operator n̂
    """
    Nvals    = np.arange(-Nmax, Nmax + 1)   # Creates a 1D array: [-Nmax, -Nmax+1, ...,-2, -1, 0, 1, 2, ..., Nmax]
    dim      = Nvals.size
    basis    = np.eye(dim)                  # Creates a (2*Nmax + 1)x(2*Nmax + 1) identity matrix

    n_charge = np.diag(Nvals.astype(float)) # Places the values of Nvals on the diagonal

    proj_terms = []
    for i, N in enumerate(Nvals):
        weighted = (N - ng)**2 * np.outer(basis[i], basis[i]) # Outer product = |N⟩⟨N| 
        proj_terms.append(weighted)
    Hc = 4 * Ec * np.sum(proj_terms, axis=0)

    hop_terms = []
    for i in range(dim - 1):
        hop_terms.append(
            np.outer(basis[i], basis[i + 1]) +
            np.outer(basis[i + 1], basis[i])
        )
    Hj = -(Ej / 2) * np.sum(hop_terms, axis=0)

    H = Hc + Hj
    return (H, n_charge)


# ============================================================
# Qubit Eigenproblem and Charge Matrix Elements
# ============================================================
def _extract_Ej_levels_and_n_phi(Lj, Csigma, Nmax):
    """
    Diagonalize the CPB Hamiltonian and extract qubit eigenenergies
    and charge matrix elements in the energy eigenbasis.

    Procedure:
    1. Construct CPB Hamiltonian in the charge basis
    2. Diagonalize Ĥ |φ_j⟩ = E_j |φ_j⟩
    3. Transform n̂ → n̂_φ = U† n̂ U

    Args:
        Lj     : Josephson inductance (H)
        Csigma : Total qubit capacitance ΣC
        Nmax   : Charge basis truncation

    Returns:
        Ej_levels : Qubit eigenenergies E_j (J)
        n_phi     : Charge operator ⟨φ_i|n̂|φ_j⟩
    """
    phi0            = hbar / (2 * e)
    Ec              = (e**2) / (2 * Csigma)
    Ej_levels       = phi0**2 / Lj
    ng              = 0.0

    H, n_charge     = _cpb_hamiltonian(Ec, Ej_levels, ng, Nmax)
    Ej_levels, U    = np.linalg.eigh(H)             # Returns eigenvalues (Ej_levels) and eigenvectors (columns of U)

    n_phi       = U.conj().T @ n_charge @ U     # Use @ for matrix multiplication

    return (Ej_levels, n_phi)

# ============================================================
# Qubit-Cavity Coupling Strength
# ============================================================
def _compute_gij(i, j, beta, Vzpf, n_phi):
    """
    Compute the qubit-cavity coupling rate g_ij between levels i and j.

    The coupling strength is:
        g_ij = (2 β e V_zpf / ℏ) ⟨φ_i|n̂|φ_j⟩

    Args:
        i     : Initial qubit level index
        j     : Final qubit level index
        beta  : Capacitive coupling ratio β
        Vzpf  : Zero-point RMS voltage V_zpf
        n_phi : Charge operator in eigenbasis

    Returns:
        g_ij : Coupling strength (rad/s)
    """
    gij = (2 * beta * e * Vzpf / hbar) * n_phi[i, j]
    return gij


# ============================================================
# Full System Hamiltonian (Local Oscillator Method)
# ============================================================
def compute_hamiltonian_lom(Lj, Lr, Cj, Cs, Cg, Cr, Nc, Nq, Nmax):
    """
    Construct the full Hamiltonian for a CPB qubit coupled to a resonator.

    The total Hamiltonian is:
        Ĥ = Ĥ_r + Ĥ_q + Ĥ_i

    with:
        Ĥ_r = ℏ ω_r ( Î_q ⊗ â† â )
  
                Nq-1
        Ĥ_q = ℏ  Σ   ω_j |φ_j⟩⟨φ_j| ⊗ Î_c
                j=0

               Nq-1
        Ĥ_i = ℏ  Σ   g_ij |φ_i⟩⟨φ_j| ⊗ ( â + â† )
              i,j=0
              
    The Hilbert space ordering is:
        |qubit_state⟩ ⊗ |cavity_state⟩

    Args:
        Lj   : Josephson inductance (H)
        Lr   : Resonator inductance (H)
        Cj   : Junction capacitance (F)
        Cs   : Shunt capacitance (F)
        Cg   : Coupling capacitance (F)
        Cr   : Resonator capacitance (F)
        Nc   : Cavity Fock space dimension
        Nq   : Number of qubit energy levels
        Nmax : Charge basis truncation

    Returns:
        Dictionary containing:
            H  : Total Hamiltonian
            Hr : Resonator Hamiltonian
            Hq : Qubit Hamiltonian
            Hi : Interaction Hamiltonian
    """
    (Csigma, beta, Wr, Vzpf, Ic_hat, Iq_hat, a_hat, a_hat_dagger) = _compute_constants(
        Lr, Cj, Cs, Cg, Cr, Nc, Nq
    )

    (Ej_levels, n_phi) = _extract_Ej_levels_and_n_phi(Lj, Csigma, Nmax)  # E = ℏ ω, therefore Ej = ℏ ω_j

    Hr          = hbar * Wr * np.kron(Iq_hat, a_hat_dagger @ a_hat)      # np.kron(x, y) = x ⊗ y

    proj_terms = []
    for j in range(Nq):
        x = Ej_levels[j] * np.kron(
            np.outer(Iq_hat[j], Iq_hat[j]),
            Ic_hat
        )
        proj_terms.append(x)
    Hq          = np.sum(proj_terms, axis=0)

    proj_terms2 = []
    for i in range(Nq):
        for j in range(Nq):
            gij = _compute_gij(i, j, beta, Vzpf, n_phi)
            x   = gij * np.kron(
                np.outer(Iq_hat[i], Iq_hat[j]),
                a_hat + a_hat_dagger
            )
            proj_terms2.append(x)
    Hi          = hbar * np.sum(proj_terms2, axis=0)

    H           = Hr + Hq + Hi
    return dict(H=H, Hr=Hr, Hq=Hq, Hi=Hi)
