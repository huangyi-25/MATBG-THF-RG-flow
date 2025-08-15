# MATBG-THF-RG-flow

## Overview

This repository contains Mathematica notebooks for studying **Magic Angle Twisted Bilayer Graphene (MATBG)** using **Topological Heavy Fermion (THF)** theory and **Renormalization Group (RG)** flow analysis. The work focuses on understanding the electronic properties and phase transitions in twisted bilayer graphene at magic angles through theoretical calculations and numerical simulations.

## Repository Contents

### 📓 Notebooks

1. **`RG_numerical.nb`** (775KB, 15,641 lines)
   - Main numerical renormalization group calculations
   - Comprehensive RG flow analysis with multiple coupling constants
   - Detailed plotting and visualization of RG flows
   - Parameter studies for various physical regimes

2. **`renormalization_chiral_limit.nb`** (349KB, 10,517 lines)
   - Specialized calculations in the chiral limit
   - Band structure analysis with 6-band basis (c1, c2, c3, c4, f1, f2)
   - Eigensystem calculations and Green's function analysis
   - Projected Green function computations

### 🐍 Python Scripts

3. **`dispersion_fcomp.py`** (5.6KB, 161 lines)
   - Band dispersion plotting with f-component visualization
   - 6-band effective Hamiltonian diagonalization
   - High-resolution momentum-space band structure calculation
   - Publication-quality matplotlib plots with Times New Roman fonts

4. **`50_dispersion_fcomp.py`** (5.4KB, 158 lines)
   - Alternative parameter set for different physical regime
   - Non-zero mass term and different coupling parameters
   - Comparative study capabilities

### 📄 Documentation

- **`LICENSE`** - MIT License for open-source distribution

## Physics Background

### Magic Angle Twisted Bilayer Graphene (MATBG)
MATBG is formed by stacking two graphene layers with a small twist angle (~1.05°). At this "magic angle," the system exhibits:
- Flat bands near the Fermi level
- Strong electronic correlations
- Superconducting and insulating phases
- Rich phase diagram with competing orders

### Key Physical Parameters

The notebooks define several important physical constants and coupling parameters:

#### Material Parameters
- **Fermi velocity** (`vs`): c-fermion Dirac velocity
- **Interlayer coupling** (`gamma`): energy separation between remote bands
- **Moire unit cell area** (`Omega0`): ~156 nm²
- **Cutoff momentum** (`Lambda0`): Wave vector at 1st mBZ boundary
- **Energy scale** (`E0`): Characteristic energy combining velocity and mass gap

#### Coupling Constants (all in meV)
- **V0** (~48.3): Density-density interaction
- **U0** (57.95): On-site Hubbard interaction  
- **W10** (44.03): Inter-sublattice coupling
- **W30** (50.20): Extended-range interaction
- **J0** (16.38): Exchange coupling
- **Jp0** (16.38): Pair-hopping coupling

#### Python Script Parameters
Two different parameter regimes are implemented:

**`dispersion_fcomp.py`** (Chiral limit):
- Mass term: M = 0.0 meV
- Pair hopping: v'* = 0.0 eV⋅Å  
- Focus on massless Dirac physics

**`50_dispersion_fcomp.py`** (Massive regime):
- Mass term: M = 3.697 meV
- Pair hopping: v'* = 1.623 eV⋅Å
- Studies gap opening and massive fermion behavior

## Theoretical Framework

### Renormalization Group Analysis
The notebooks implement a comprehensive RG flow study that:

1. **Tracks coupling evolution** from high-energy cutoff (`Lambda0`) to low energies
2. **Identifies fixed points** and phase transitions
3. **Analyzes stability** of different electronic phases
4. **Studies crossover behavior** between different regimes

### Band Structure Calculation
The chiral limit notebook includes:
- **6-band effective model** for the moire superlattice
- **Eigensystem analysis** of the continuum Hamiltonian
- **Momentum-resolved** band structure calculations
- **Green's function** formalism for many-body effects

## Usage

### Prerequisites

#### For Mathematica Notebooks:
- **Mathematica 12.3** or later
- **MaTeX package** (for LaTeX-quality plots)

#### For Python Scripts:
- **Python 3.7+**
- **NumPy** (`pip install numpy`)
- **Matplotlib** (`pip install matplotlib`)

### Running the Code

#### 1. Clone the Repository
```bash
git clone https://github.com/huangyi-25/MATBG-THF-RG-flow.git
cd MATBG-THF-RG-flow
```

#### 2. Mathematica Notebooks
- Launch Mathematica
- Open either `RG_numerical.nb` or `renormalization_chiral_limit.nb`
- Execute cells sequentially or use "Evaluate Notebook"

#### 3. Python Band Structure Plotting
```bash
# Install dependencies
pip install numpy matplotlib

# Run the main dispersion script
python dispersion_fcomp.py

# Or run the alternative parameter set
python 50_dispersion_fcomp.py
```

#### 4. Parameter Modification
- **Mathematica**: Physical parameters are defined in the initial cells
- **Python**: Edit the constants section at the top of the script:
  ```python
  vstar = -4.303  # eV * angstrom (Fermi velocity)
  gamma = -24.75  # meV (interlayer coupling)
  M = 0.0         # meV (mass term)
  theta0 = 2 * np.pi * 1.05 / 360  # Magic angle
  ```

### Key Computational Features

#### RG Flow Calculations (`RG_numerical.nb`)
- **Numerical integration** of RG equations
- **Multi-parameter tracking** (8 coupling constants)
- **Flow visualization** with high-quality plots
- **Critical point analysis**

#### Band Structure Analysis (`renormalization_chiral_limit.nb`)
- **Analytical diagonalization** of effective Hamiltonian
- **Momentum-space** eigenvalue calculations
- **Symmetry analysis** in different geometric configurations
- **Green's function** spectral properties

#### Python Band Dispersion (`dispersion_fcomp.py`)
- **6-band effective model** implementation using NumPy
- **Momentum path**: Γ_M → K_M → M_M → Γ_M through moire Brillouin zone
- **f-component analysis** showing hybridization with heavy fermions
- **High-resolution plotting** (80 k-points along each segment)
- **Physical parameters**:
  - Magic angle: θ₀ = 1.05°
  - Moire lattice constant: a_M ≈ 100 nm
  - Energy scale: ~150 meV range
- **Output**: Publication-ready PDF plots with f-component weight visualization

## Results and Applications

The notebooks enable investigation of:

- **Phase diagrams** as functions of twist angle and filling
- **Correlation effects** in the flat band regime  
- **Superconducting instabilities** and pairing mechanisms
- **Topological properties** and band topology
- **Transport signatures** and experimental observables

## Citation

If you use this code in your research, please cite:

```bibtex
@software{matbg_thf_rg_2025,
  author = {huangyi-25},
  title = {MATBG-THF-RG-flow: Renormalization Group Analysis of Magic Angle Twisted Bilayer Graphene},
  year = {2025},
  url = {https://github.com/huangyi-25/MATBG-THF-RG-flow}
}
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add your modifications
4. Submit a pull request with detailed description

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaborations, please open an issue in the GitHub repository.

---

*This work contributes to the understanding of emergent phenomena in moiré quantum materials and correlated electron systems.* 
