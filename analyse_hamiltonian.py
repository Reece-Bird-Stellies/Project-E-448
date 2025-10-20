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
    # Protection: check if H is valid
    H_0 = {
            'qubit_frequency_ghz': 0,
            'cavity_frequency_ghz': 0,
            'anharmonicity_mhz': 0,
            'coupling_g_mhz': 0,
            'kappa_khz': 0,
            'h_dimension': 0
        }
    if H is None:
        print("Hamiltonian is None - returning zero results")
        return H_0
    
    if not isinstance(H, np.ndarray):
        print(f"Hamiltonian is not an array (type: {type(H)}) - returning zero results")
        return H_0
    
    if not hasattr(H, 'shape') or len(H.shape) != 2:
        print(f"Hamiltonian is not a 2D array - returning zero results")
        return H_0
    
    if H.shape[0] < 3 or H.shape[1] < 3:
        print(f"Hamiltonian is too small ({H.shape}) - need at least 3x3 - returning zero results")
        return H_0
    
    try:
        E_matrix = _test_and_diagonalize(H)

        F0 = E_matrix[0, 0] / h
        F1 = E_matrix[1, 1] / h
        F2 = E_matrix[2, 2] / h

        f12 = F2 - F1
        f01 = F1 - F0
        f20 = F2 - F0

        qubit_frequency     = (f01) * 1e-9  # in GHz
        cavity_frequency    = (f20) * 1e-9  # in GHz
        ana_harm            = (f12 - f01) * 1e-6  # in MHz
        
        g = 0
        kappa = 0
        return {
            'qubit_frequency_ghz': qubit_frequency,
            'cavity_frequency_ghz': cavity_frequency,
            'anharmonicity_mhz': ana_harm,
            'coupling_g_mhz': g,
            'kappa_khz': kappa,
            'h_dimension': H.shape[0]
        }
    
    except Exception as e:
        print(f"Error analyzing Hamiltonian: {e}")
        return H_0
