[![Checks](https://github.com/aladshaw3/porous_flow_matlab/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/aladshaw3/porous_flow_matlab/actions/workflows/unit_tests.yml)
[![codecov](https://codecov.io/gh/aladshaw3/porous_flow_matlab/branch/main/graph/badge.svg)](https://codecov.io/gh/aladshaw3/porous_flow_matlab) 

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aladshaw3/porous_flow_matlab)

## WIP: Use with Caution ##


# Porous Flow in MATLAB
This repository provides a MATLAB interface for solving 2D PDEs resulting 
from reactive water flow through a porous media. The physics include:

 - Pressure driven advective/convective flow

 - Thermal balance of the fluid phase 

 - Mass balance for chemicals in the fluid phase

# Requirements

 - MATLAB
 - [PDE Toolbox](https://www.mathworks.com/help/pde/)


# Equation Set

 - Pressure driven Darcy flow

$$ \nabla{ ( K  \cdot \nabla{P} ) } = 0 $$

$$ \vec{v} = -K \cdot \nabla{P} $$

 - Thermal energy balance

$$ (\varepsilon \rho c_{pw} + (1- \varepsilon) \rho_{s} c_{ps}) \frac{\partial T}{\partial t} - \nabla{ ( (\varepsilon K_w + (1- \varepsilon) K_s ) \cdot \nabla{T} ) } = - (\varepsilon \rho c_{pw}) \vec{v} \cdot \nabla{T} + \varepsilon \sum{\Delta H_j \cdot r_j } $$

 - Mass balances (for each species)

$$ \varepsilon \frac{\partial C_i}{\partial t} - \nabla{ ( \delta_i \varepsilon D_i \cdot \nabla{C_i} ) } = - \delta_i \varepsilon \vec{v} \cdot \nabla{C_i} + \varepsilon \sum{u_j \cdot r_j }$$

 - Reaction term

$$ r_j = k_j \cdot exp( - \frac{E_j}{R T} ) \cdot \prod{C_l^{s_l}} $$

 - Mobility term

If species is mobile...

$$ \delta_i = 1 $$

Otherwise...

$$ \delta_i = 0 $$

# Citation 

Ladshaw, A.P., "Porous Flow MATLAB: A simple 2D interface for modeling flow and reactions
through porous media," https://github.com/aladshaw3/porous_flow_matlab, 
Accessed (Month) (Day), (Year).
