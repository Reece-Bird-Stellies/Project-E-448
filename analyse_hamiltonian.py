import numpy as np
from scipy.constants import Planck, e, hbar, h

def _test_and_diagonalize(H):
    if np.allclose(H, H.conj().T):
        print("Hamiltonian is Hermitian")
        evals, evecs    = np.linalg.eigh(H)             # diagonalize if you want energies/eigenstates
        E_matrix        = np.diag(evals)                # diagonalized Hamiltonian

        # Verify diagonalization: V^† H V = Λ
        H_diag          = evecs.conj().T @ H @ evecs
        print("Is diagonal?", np.allclose(H_diag, E_matrix))

        # (Optional) reconstruct original: H = V Λ V^†
        H_back          = evecs @ E_matrix @ evecs.conj().T
        print("Reconstruction ok?", np.allclose(H_back, H))
    else:
        print("Hamiltonian is NOT Hermitian")
    return E_matrix


def analyse_hamilotonian(Hamiltonians, Nc, Nq, analyse_type="LOM"):
    """Analyse properties of Hamiltonian H."""
    # Protection: check if H is valid
    # qubit first ordering
    try: 
        H  = Hamiltonians['H']  # Extract the Hamiltonians from the dictionary
        Hq = Hamiltonians['Hq']
        Hr = Hamiltonians['Hr']
    except:
        H = Hamiltonians
    
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
        if analyse_type == "EPR":
            E_matrix   = _test_and_diagonalize(H)
            
            F = []
            for energy in np.diag(E_matrix):
                F.append((energy / h)*1e-9)   # Convert E matrix to frequencies in GHz and store in frequency list
                
            f21 = F[2] - F[1]
            f10 = F[1] - F[0]
            f20 = F[2] - F[0]
            
            g                   = abs((H[0, 2*Nc] / h) * 1e-6)  # in MHz
            qubit_frequency     = f10
            cavity_frequency    = f20
            ana_harm            = (f21 - f10) * 1e3  # in MHz
        else:
            E_matrix   = _test_and_diagonalize(H)
            E_q_matrix = _test_and_diagonalize(Hq)
            E_r_matrix = _test_and_diagonalize(Hr)
            
            F_q = []
            for energy in np.diag(E_q_matrix):
                F_q.append((energy / h)*1e-9)   # Convert E matrix to frequencies in GHz and store in frequency list
            
            F_r = []
            for energy in np.diag(E_r_matrix):
                F_r.append((energy / h)*1e-9)   # Convert E matrix to frequencies in GHz and store in frequency list
                

            qubit_frequency     = F_q[Nc] - F_q[0]
            cavity_frequency    = F_r[Nq] - F_r[0]
            ana_harm            = ((F_q[2*Nc] - F_q[Nc]) - (F_q[Nc] - F_q[0])) *1e3
            g                   = abs((H[1, Nc] / h) * 1e-6)
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
