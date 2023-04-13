[![Checks](https://github.com/aladshaw3/porous_flow_matlab/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/aladshaw3/porous_flow_matlab/actions/workflows/unit_tests.yml)
[![codecov](https://codecov.io/gh/aladshaw3/porous_flow_matlab/branch/main/graph/badge.svg)](https://codecov.io/gh/aladshaw3/porous_flow_matlab) 

# Porous Flow in MATLAB
This repository provides a MATLAB interface for solving 2D PDEs resulting 
from reactive water flow through a porous media. The physics include:

 - Pressure driven advective/convective flow

 - Thermal balance of the fluid phase 

 - Mass balance for chemicals in the fluid phase


# Equation Set

 - Pressure driven Darcy flow

$$ \nabla{ ( K  \cdot \nabla{P} ) } = 0 $$

$$ \vec{v} = -K \cdot \nabla{P} $$

 - Thermal energy balance

$$ (\varepsilon \rho c_{pw} + (1- \varepsilon) \rho_{s} c_{ps}) \frac{\partial T}{\partial t} - \nabla{ ( (\varepsilon K_w + (1- \varepsilon) K_s ) \cdot \nabla{T} ) } = - (\varepsilon \rho c_{pw}) \vec{v} \cdot \nabla{T} + \varepsilon \sum{\Delta H_j \cdot r_j } $$

 - Mass balances (for each species)

$$ \varepsilon \frac{\partial C_i}{\partial t} - \nabla{ ( D_i \cdot \nabla{C_i} ) } = - \varepsilon \vec{v} \cdot \nabla{C_i} + \varepsilon \sum{u_j \cdot r_j }$$

 - Reaction term

$$ r_j = k_j \cdot exp( - \frac{E_j}{R T} ) \cdot \prod{C_i^{s_i}} $$