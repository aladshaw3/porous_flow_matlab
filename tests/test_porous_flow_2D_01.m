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

% Calculation of the effective thermal conductivity of water [W/m/k]
Kew = obj.EffectiveThermalConductivityWater(101350, 298, 0.1);
assert( abs(Kew-9.36121272563573) < 1e-6)

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


%% Test 2 - Basic model setup

% Create geometry for the pde
%   2D geometry 

% Rectangle is code 3, 4 sides (in m)
R1 = [3,4, 0,4, 4,0, 0,0, 4,4]';
geom = R1;

% Create geometry
% E1 - bottom, E2 - right, E3 - top, E4 - left
g = decsg(geom);

obj = porous_flow_2D();

obj.set_geometry_from_edges(g);

dmat = @(location,state) obj.d_coeff_fun(1,location,state);
cmat = @(location,state) obj.c_coeff_fun(1,location,state);

% Specify what the model coefficients are
specifyCoefficients(obj.model,"m",0,"d",dmat,"c",cmat,"a",0,"f",[0;0;0]);

% BC Formats
%
%       Dirichlet:  h*u=r
%
%       Neumann:    n * (c * grad(u)) = g - q*u

% Apply BCs
% Enter
he = [0 0 0; 
     0 1 0; 
     0 0 1];

re = [0;300;1];

qe = [0 0 0; 
     0 0 0; 
     0 0 0];

ge = [1;0;0];

applyBoundaryCondition(obj.model,"mixed", ...
                             "Edge",3, ...
                             "h",he,"r",re,"q",qe,"g",ge);


% exit
hw = [1 0 0; 
     0 0 0; 
     0 0 0];

rw = [101350;0;0];

qw = [0 0 0; 
      0 0 0; 
      0 0 0];

gw =  [0; 0; 0];

applyBoundaryCondition(obj.model,"mixed", ...
                             "Edge",1, ...
                             "h",hw,"r",rw,"q",qw,"g",gw);


% Set initial conditions 
u0 = [101350;298;0];
setInitialConditions(obj.model,u0);
generateMesh(obj.model);

% Solve the model
obj.model.SolverOptions.ReportStatistics = 'on';
obj.model.SolverOptions.AbsoluteTolerance = 1e-4; % ODE opt
obj.model.SolverOptions.RelativeTolerance = 1e-4; % ODE opt
obj.model.SolverOptions.ResidualTolerance = 1e-6; % Nonlinear opt
obj.model.SolverOptions.MaxIterations = 30;       % Nonlinear opt
obj.model.SolverOptions.MinStep = 0.001;          % Min step size 
obj.model.SolverOptions.ResidualNorm = 2;         % L-2 norm
obj.model.SolverOptions.MaxShift = 500;           % Lanczos solver shift
obj.model.SolverOptions.BlockSize = 50;           % Block size for Lanczos recurrence

nsteps=10;
t_span = linspace(0,100,nsteps);
results = solvepde(obj.model,t_span);
u = results.NodalSolution;

f1 = figure;
pdeplot(obj.model,"XYData", u(:,1,end),"ZData",u(:,1,end) ,Mesh="on", ColorMap="jet")
% saveas(f1,'output/Gifs/pressure_test01.png');
close(f1);

% Create gif for temperature
first_command = 'pdeplot(model,"XYData",u(:,2,INDEX),"ZData",u(:,2,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,2,:)));','umin = min(min(u(:,2,:)));',...
                'axis([0 4 0 4 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel T;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output2 = create_gif(command_set,variable_set,index_limit,1,'temp_solution_test01');


% Create gif for concentration 
first_command = 'pdeplot(model,"XYData",u(:,3,INDEX),"ZData",u(:,3,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,3,:)));','umin = min(min(u(:,3,:)));',...
                'axis([0 4 0 4 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel C;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output3 = create_gif(command_set,variable_set,index_limit,1,'conc_solution_test01');
