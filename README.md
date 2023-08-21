[![Checks](https://github.com/aladshaw3/porous_flow_matlab/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/aladshaw3/porous_flow_matlab/actions/workflows/unit_tests.yml)
[![codecov](https://codecov.io/gh/aladshaw3/porous_flow_matlab/branch/main/graph/badge.svg)](https://codecov.io/gh/aladshaw3/porous_flow_matlab) 

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aladshaw3/porous_flow_matlab)

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

> If species is mobile...

$$ \delta_i = 1 $$

> Otherwise...

$$ \delta_i = 0 $$

# Examples

All examples discussed have their source code under the `tests` directory.

 ## (1) Basic Usage

Test file `tests/test_porous_flow_2D_01.m` contains some basic calculation 
checks for material properties, but also exercises the simulation tools basic 
use cases. The test case involves a simple square domain (4m x 4m) and a 
single chemical species with 2 first-order loss reactions. 

### Creating an instance
 
This example starts by initiallizing the `porous_flow_2D` object with vector
spaces for 1 species, 2 reactions, and a single subdomain.

```
obj = porous_flow_2D(1, 2, 1);
```

### Importing geometry for meshing

The geometry is created using the built in Matlab `decsg` function, however, 
you can import any geometry file to use. Once you have a defined geometry, 
you pass it into the `porous_flow_2D` object. You can also define the 
`GeometryOrder` and `Hmax` (i.e., maximum element size). The Matlab PDE toolbox
will automatically generate the mesh from this information. 

```
obj.set_geometry_from_edges(g, "quadratic", 0.25);
```

### Editable material properties

Next, we can setup the information and parameters for the chemical species 
and reactions (as well as other properties such as particle size, bulk densities,
thermal properties, etc.). A full listing of these properties can be found 
by reading through the `porous_flow_2D` object itself, but are also listed 
below. 

```
    %% Public properties
    properties (Access = public)
        particle_diameter     % particle diameter in m
        particle_porosity     % Volume of micro-voids per bulk solids volume
        muA                   % Viscosity pre-exponential in kg/m/s
        muB                   % Viscosity temperature relation B in K
        muC                   % Viscosity temperature relation C in K 
        rhoA                  % Density temperature relation A in K^-2
        rhoB                  % Density temperature relation B in K^-1
        rhoC                  % Density temperature relation C in -
        rhoRef                % Reference Density in kg/m^3
        bulk_porosity         % system bulk porosity (vol voids / vol bulk solids)
        KozenyCarmannConst    % Kozeny-Carmann coefficient (dimensionless)
        diff_factor           % Diffusivity factor 
        alpha                 % Dispersion length factor in m
        refDiff               % Reference diffusivity in m^2/s
        refDiffTemp           % Reference diffusivity temperature in K
        cpA                   % Specific heat temperature parameter in K^-2
        cpB                   % Specific heat temperature parameter in K^-1
        cpC                   % Specific heat reference parameter -
        KwA                   % Thermal conductivity of water slope in K^-1
        KwB                   % Thermal conductivity of water intercept -
        char_len              % A characteristic length for the domain in m
        bulk_solids_dens      % Bulk solids density in the domain in kg/m^3
        bulk_solids_Cp        % Bulk solids heat capacity in J/kg/K
        bulk_solids_cond      % Bulk solids thermal conductivity in J/m/K/s

        mobile_spec_idx       % Matrix of 1s and 0s (where 1 indicates the species is mobile)
        rxn_stoich            % Stoichiometry matrix
        rxn_powers            % Matrix for reaction power terms
        rxn_rate_const        % Matrix for reaction rate constants (units will vary)
        rxn_act_energy        % Matrix for reaction rate activation energies (J/mol)
        rxn_enthalpy          % Matrix for reaction enthalpies (J/mol)

        model                 % PDE model object 

        boundaries            % Object to identify boundaries 
    end
```

### Setting up FEM coefficient matrices

This Matlab object utilizes the [General PDE](https://www.mathworks.com/help/pde/pde-problem-setup.html) 
toolbox system in Matlab. This system uses a faily complex, but comprehensive way
to setup the matrix system for all variables in your system simultaneously. More 
details can be found in the links below:

 - [c matrix coefficients](https://www.mathworks.com/help/pde/ug/c-coefficient-for-systems-for-specifycoefficients.html)
 - [f vector coefficients](https://www.mathworks.com/help/pde/ug/f-coefficient-for-specifycoefficients.html)
 - [m, d, and a coefficients](https://www.mathworks.com/help/pde/ug/m-d-or-a-coefficient-for-systems.html)

To simplify this process, after you have finished editing any of the material 
properties from `porous_flow_2D` object, you can simply call the `set_coefficients`
function to automatically have all these coefficient matrices built automatically.

**NOTE**: You will want to do this BEFORE applying Boundary Conditions (BCs)

```
obj.set_coefficients();
```

### Setting up Boundary Conditions

Just like the coefficient matrices, the Matlab toolbox has a relatively complex
set of subroutines for [applying BCs](https://www.mathworks.com/help/pde/ug/steps-to-specify-a-boundary-conditions-object.html).
To simplify this, the `porous_flow_2D` object will automatically apply all the 
proper BCs for your problem based on whether it is an "input" or "output" boundary.

#### Input Boundary

To set the input boundaries, you need to supply:

 - The set of boundary ids
 - A fluid velocity at each id
 - An input temperature at each id
 - An input concentration for each species at each id (NxID)

```
inbound_set = [3,4];
velocity_set = [0.005,0.0075];
temperature_set = [298,298];
concentration_matrix = [1,0.5];
obj.set_input_boundaries(inbound_set, velocity_set, ...
                    temperature_set, concentration_matrix);
```

#### Output Boundary

To set the output boundaries, you only need to supply:

 - The set of boundary ids
 - Reference state pressure at the output boundaries 

```
outbound_set = [1];
pressure_set = [101350];
obj.set_output_boundaries(outbound_set, pressure_set);
```

**NOTE**: Any boundary ID where you DO NOT specify boundary conditions will 
utilize the so-called "Natural BC" for a FEM problem, which is equivalent 
to a Neumann-type BC where the slope of the variable at the bounary is zero.


### Setting up Initial Conditions

To set initial conditions (ICs), you need to provide values for each of the variables
within each subdomain of the solution space. In this example, we only have 1 
subdomain, so we just provide the variable ICs for pressure, temperature, and 
the 1 chemical species.

```
subdomain_set = [1];
pressure_set = [101350];
temperature_set = [298];
concentration_matrix = [0];
obj.set_initial_conditions(subdomain_set,pressure_set,temperature_set, concentration_matrix);
```

### Setting up the time span and solving

To solve the system, simply create a `t_span` for the number of time steps 
and time positions to solve for, then call the `solve_system` function. At output,
you will recieve "results" as a [TimeDependentResults](https://www.mathworks.com/help/pde/ug/pde.timedependentresults.html) 
object. From this object, you can extract things such as the "NodalSolution" for 
post-processing and generating plots. 

```
nsteps=20;
t_span = linspace(0,200,nsteps);
results = obj.solve_system(t_span);
u = results.NodalSolution;
```

### Example 1 - Results

![Fig1](output/Gifs/pressure_test01.png?raw=true)

**Figure 1**: Pressure distribution (Pa)

![Fig2](output/Gifs/temp_solution_test01.gif?raw=true)

**Figure 2**: Temperature dynamics (K)


![Fig3](output/Gifs/conc_solution_test01.gif?raw=true)

**Figure 3**: Concentration dynamics (mol/m^3)

---

 ## (2) Injection well and Draw well

---

 ## (3) Multi-domain Properties and Chemical Breakthrough

---

# Citation 

Ladshaw, A.P., "Porous Flow MATLAB: A simple 2D interface for modeling flow and reactions
through porous media," https://github.com/aladshaw3/porous_flow_matlab, 
Accessed (Month) (Day), (Year).
