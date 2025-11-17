# 🌞 Plasma Convection Simulation

2D hydrodynamic simulation of convective motions in stellar interiors, developed at the **Institute of Theoretical Astrophysics, University of Oslo**.

## 📌 Features
- **Custom hydrodynamic solver** implementing continuity, momentum and energy equations
- **Finite difference methods** with upwind scheme for numerical stability
- **Multi-parameter analysis**: temperature, velocity, pressure, and density evolution
- **Variable adiabatic gradient (∇)** to study convective instability thresholds
- **Dynamic time stepping** with CFL condition for simulation stability
- **Gaussian perturbations** to trigger convective flows in hydrostatic equilibrium

## 🔧 Technical Implementation

### Custom Code
- **`hydro_solver.py`** - Main hydrodynamic solver (original implementation)
- Boundary condition handling and initial state setup
- Numerical schemes for fluid equations

### Visualization Framework
- **`FVis3.py`** - Visualization module from University of Oslo
- Used for animation and data visualization
- **FVis3 Repository**: [lars-frogner/FVis on GitHub](https://github.com/lars-frogner/FVis)

## 🎬 Simulation Results

### Temperature Evolution
| ∇ = 0.4 (Stable) | ∇ = 2.0 (Convective) | ∇ = 5.0 (Turbulent) |
|------------------|---------------------|---------------------|
| [∇=0.4.zip](simulations/temperature/vid_T-nabla-0.4.mp4.zip) | [∇=2.0.zip](simulations/temperature/vid_T-nabla-2.mp4.zip) | [∇=5.0.zip](simulations/temperature/vid_T-nabla-5.mp4.zip) |

### Vertical Velocity Patterns
| ∇ = 0.4 | ∇ = 2.0 | ∇ = 5.0 |
|---------|---------|---------|
| [V0.4.zip](simulations/vertical_velocity/vid_w-nabla-0.4.mp4.zip) | [V2.0.zip](simulations/vertical_velocity/vid_w-nabla-2.mp4.zip) | [V5.0.zip](simulations/vertical_velocity/vid_w-nabla-5.mp4.zip) |

### Pressure Distribution
| ∇ = 0.4 | ∇ = 2.0 | ∇ = 5.0 |
|---------|---------|---------|
| [P0.4.zip](simulations/pressure/vid_pressure-nabla-0.4.mp4.zip) | [P2.0.zip](simulations/pressure/vid_pressure-nabla-2.mp4.zip) | [P5.0.zip](simulations/pressure/vid_pressure-nabla-5.mp4.zip) |

### Density & Energy Variations
| ∇ = 0.4 | ∇ = 2.0 | ∇ = 5.0 |
|---------|---------|---------|
| [ρ0.4.zip](simulations/density/vid_energy_density-nabla-0.4.mp4) | [ρ2.0.zip](simulations/density/vid_energy_density-nabla-2.mp4.zip) | [ρ5.0.zip](simulations/density/vid_energy_density-nabla-5.mp4.zip) |

## 📊 Key Findings
- **∇ < 0.4**: System remains in hydrostatic equilibrium, no convection develops
- **∇ > 0.4**: Convective cells emerge naturally with well-defined patterns
- **∇ = 2.0**: Optimal convection with clear cellular structures
- **∇ = 5.0**: Strong instabilities lead to turbulent behavior and numerical noise
- **Energy conservation**: System maintains physical consistency across all parameters

## 🛠️ Installation & Usage

### Prerequisites
- Python 3.7+
- NumPy
- Matplotlib

### Setup
```bash
# Clone repository
git clone https://github.com/Ignacioferris/plasma-convection-simulation.git
cd plasma-convection-simulation

# Run main simulation
python src/hydro_solver.py
