import numpy as np
import matplotlib.pyplot as plt

# Set font to Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'

# Define constants (units indicated in comments)
vstar = -4.303  # eV * angstrom
# vstarp = 1.623  # eV * angstrom
vstarp = 0.0  # eV * angstrom
# M = 3.697  # meV
M = 0.0  # meV
gamma = -24.75  # meV

theta0 = 2 * np.pi * 1.05 / 360
KDirac = 1.703 * 1000  # 1/(1000 angstrom) = 1/(100 nm)

# Length of wave vectors q_j
ktheta = 2 * KDirac * np.sin(theta0 / 2)  # 1/(1000 angstrom) = 1/(100 nm)

# Moiré lattice constant
aM = (4 * np.pi) / (3 * ktheta)  # 1000 angstrom = 100 nm

# Characteristic length scale
lambda_ = 0.3375 * aM  # 1000 angstrom = 100 nm

# Define Pauli matrices
zero0 = np.array([[0, 0],
                  [0, 0]])

sigma0 = np.array([[1, 0],
                   [0, 1]])

sigmax = np.array([[0, 1],
                   [1, 0]])

sigmay = np.array([[0, -1j],
                   [1j, 0]])

sigmaz = np.array([[1, 0],
                   [0, -1]])

def hc(kx, ky, eta, vstar, M):
    """Single-particle Hamiltonian for the c fermions"""
    block11 = zero0
    block12 = vstar * (eta * kx * sigma0 + 1j * ky * sigmaz)
    block21 = vstar * (eta * kx * sigma0 - 1j * ky * sigmaz)
    block22 = M * sigmax
    return np.block([[block11, block12],
                     [block21, block22]])

def hcf(kx, ky, eta, gamma, vstarp, lambda_):
    """Coupling between c and f fermions"""
    prefactor = np.exp(-((kx**2 + ky**2) * lambda_**2) / 2)
    block11 = gamma * sigma0 + vstarp * (eta * kx * sigmax + ky * sigmay)
    block12 = zero0
    return prefactor * np.block([[block11, block12]])

def h0(kx, ky, eta):
    """Full single-particle Hamiltonian that couples f and c fermions"""
    h_c = hc(kx, ky, eta, vstar, M)
    h_cf = hcf(kx, ky, eta, gamma, vstarp, lambda_)
    h_fc = h_cf.conj().T  # Take the Hermitian conjugate
    h_ff = zero0
    return np.block([[h_c, h_fc],
                     [h_cf, h_ff]])


# Parameters
Lk = 80  # The line from Gamma_M to K_M is equally divided into Lk pieces.
dk = ktheta / Lk  # The distance between two adjacent k-points.

g2 = np.array([np.sqrt(3)/2, 1/2])
g3 = np.array([np.sqrt(3)/2, -1/2])  # g2 and g3 are unit vectors pointing from Gamma_M to K_M.

# Generate all possible coordinates (m2, m3) satisfying m2+m3<=4 && m2-m3>=0.
ContourBZCoordinatesList = (
    [[m2, 0] for m2 in range(Lk, -1, -1)] +
    [[m2, m2] for m2 in range(1, Lk//2 + 1)] +
    [[Lk//2 + m, Lk//2 - m] for m in range(1, Lk//2 + 1)]
)

# Compute corresponding k-points
ContourBZkpointList = np.array(
    dk * np.array([[m2] * np.array(g2) for m2 in range(Lk, -1, -1)] +
                   [[m2] * (np.array(g2) + np.array(g3)) for m2 in range(1, Lk//2 + 1)] +
                   [(Lk//2 * (np.array(g2) + np.array(g3)) + m * (np.array(g2) - np.array(g3))) for m in range(1, Lk//2 + 1)])
)
# Compute k-axis list
kaxisList = [0] + list(np.cumsum([np.linalg.norm(ContourBZkpointList[i] - ContourBZkpointList[i-1]) for i in range(1, len(ContourBZkpointList))]))


# Compute eigenvalues and eigenvectors, then sort them
contourBandsList = []
fComponentList = []

for i in range(len(ContourBZkpointList)):
    kx, ky = ContourBZkpointList[i]
    eigenvalues, eigenvectors = np.linalg.eigh(h0(kx, ky, 1))
    
    # Sort eigenvalues and eigenvectors
    sorted_indices = np.argsort(eigenvalues)
    sortedEigenvalues = eigenvalues[sorted_indices]
    sortedEigenvectors = eigenvectors[:, sorted_indices]
    
    contourBandsList.append(sortedEigenvalues)
    fComponentList.append([sum(abs(sortedEigenvectors[4:6,j])**2) for j in range(6)])

# Compute k-axis bands list
contourBandskaxisList = [[kaxisList[i], contourBandsList[i][j]] for j in range(len(contourBandsList[0])) for i in range(len(kaxisList))]

# Compute f-component bands k-axis list
fComponentBandskaxisList = [[kaxisList[i], contourBandsList[i][j], fComponentList[i][j]] for j in range(len(contourBandsList[0])) for i in range(len(kaxisList))]

# Compute tick locations
xticks = [0,
    np.linalg.norm(Lk * g2 * dk),
          np.linalg.norm(Lk * g2 * dk) + np.linalg.norm(Lk/2 * (g2 + g3) * dk),
          np.linalg.norm(Lk * g2 * dk) + np.linalg.norm(Lk/2 * (g2 + g3) * dk) + np.linalg.norm(Lk/2 * (g2 + g3) * dk - Lk * g2 * dk)]
yticks = np.arange(-150, 151, 50)

# Plot using Mathematica-style aesthetics
plt.figure(figsize=(2.5, 2.5))  # ImageSize -> {160, 150}
#plt.axhline(y=0, color='black', linewidth=0.8)
#plt.grid(True, linestyle='--', linewidth=0.5)
# Plot scatter points and connecting lines

for x, y, size in fComponentBandskaxisList:
    plt.scatter(x, y, s=size * 10, color='r', alpha=0.5, linewidth=0.6)
    
contourBandsMatrix = np.vstack(contourBandsList)

for j in range(6):
    x_val = kaxisList
    y_val = contourBandsMatrix[:,j]
    plt.plot(x_val, y_val, color='k', linewidth=0.75)  # Line plot for connecting points

# Formatting
plt.xlim(0, xticks[-1])
plt.ylim(-150, 150)
plt.xticks(xticks, ['', '', '', ''])
plt.yticks(yticks, yticks)
plt.gca().set_aspect('auto')

# Add ticks and spines on all sides
plt.gca().spines['top'].set_visible(True)
plt.gca().spines['right'].set_visible(True)
plt.tick_params(axis='x', which='both', top=True, labeltop=False)
plt.tick_params(axis='y', which='both', right=True, labelright=False)
plt.gca().spines['top'].set_linewidth(0.75)
plt.gca().spines['right'].set_linewidth(0.75)
plt.gca().spines['bottom'].set_linewidth(0.75)
plt.gca().spines['left'].set_linewidth(0.75)

# Save the figure as a PDF
plt.savefig("fComponentBandskaxisList.pdf", format="pdf", bbox_inches="tight")
plt.show()

print(fComponentList)
# %%
