% Basic unit tests for the porous flow object

%% Test 1 - Basic calculations with default options
obj = porous_flow_2D();


% Calculation of Viscosity (mu) [kg/m/s]
mu = obj.ViscosityWater();
assert( abs(mu-8.943720611911905e-04) < 1e-6)

% Calculation of Density (rho) [kg/m^3]
rho = obj.DensityWater();
assert( abs(rho-1.046040050086146e+03) < 1e-6)

% Calculation of wall heat transfer coefficient [W/m^2/K]
hw = obj.WallHeatTransferWater();
assert( abs(hw-3.608516299511982e+02) < 1e-6)

% Calculation of the thermal conductivity of water [W/m/k]
Kw = obj.ThermalConductivityWater();
assert( abs(Kw-0.601419383251997) < 1e-6)

% Calculation of specific heat of water [J/kg/K]
cpw = obj.SpecificHeatWater();
assert( abs(cpw-4.187121392561558e+03) < 1e-6)

% Calculation of diffusion in water [m^2/s]
D = obj.DiffusionWater();
assert( abs(D-2.288292227481152e-09) < 1e-6)

% Calculation of dispersion in water [m^2/s]
Dp = obj.DispersionWater();
assert( abs(Dp-2.288292227481152e-09) < 1e-6)

% Calculation of effective dispersion in porous media [m^2/s]
Deff = obj.EffectiveDispersionWater();
assert( abs(Deff-1.618066951388393e-09) < 1e-6)

% Calculation of Kozeny-Carmann coefficient [m^3*s/kg]
K = obj.KozenyCarmannDarcyCoeffient();
assert( abs(K-0.002518249786618) < 1e-6)

% Calculation of time-coefficient
coeff = obj.TempTimeCoeff();
assert( abs(coeff-3.899448335595932e+06) < 1e-6)

