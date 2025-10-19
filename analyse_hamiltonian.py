import numpy as np
from scipy.constants import Planck, e, hbar, h

def _test_and_diagonalize(H):
    if np.allclose(H, H.conj().T):
        print("Hamiltonian is Hermitian")
        evals, evecs    = np.linalg.eigh(H)   # diagonalize if you want energies/eigenstates
        E_matrix        = np.diag(evals)              # diagonalized Hamiltonian

        # Verify diagonalization: V^† H V = Λ
        H_diag          = evecs.conj().T @ H @ evecs
        print("Is diagonal?", np.allclose(H_diag, E_matrix))

        # (Optional) reconstruct original: H = V Λ V^†
        H_back          = evecs @ E_matrix @ evecs.conj().T
        print("Reconstruction ok?", np.allclose(H_back, H))
    else:
        print("Hamiltonian is NOT Hermitian")
    return E_matrix

def analyse_hamilotonian(H):
    """Analyse properties of Hamiltonian H."""
    E_matrix = _test_and_diagonalize(H)

    F0 = E_matrix[0] / h
    F1 = E_matrix[1] / h
    F2 = E_matrix[2] / h
    F7 = E_matrix[7] / h
