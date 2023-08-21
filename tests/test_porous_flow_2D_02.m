%% Test 3 - Injection Well Simulation

% Rectangle is code 3, 4 sides,
% followed by x-coordinates and then y-coordinates (in m)
R1 = [3,4,-1.5*10,1.5*10,1.5*10,-1.5*10,-1*10,-1*10,1*10,1*10]';
% Circle is code 1, center (5,0), radius .5 m 
C1 = [1,5,0,.05*10]';
C2 = [1,-5,0,.05*10]';
% Pad C1 with zeros to enable concatenation with R1
C1 = [C1;zeros(length(R1)-length(C1),1)];
C2 = [C2;zeros(length(R1)-length(C2),1)];
geom = [R1,C1,C2];

% Names for the two geometric objects
ns = (char('R1','C1','C2'))';

% Set formula
sf = 'R1 - C1 - C2';

% Create geometry
g = decsg(geom,sf,ns);


% Create instance of PDE object 

%                   2 chemical species, 2 reactions, 1 subdomain
%
%       Rxn:    A --> B
%               B --> A
%               Species A is mobile, Species B is not
obj = porous_flow_2D(2, 2, 1);

% Set geometry based on edges and setup functions
obj.set_geometry_from_edges(g,"linear",0.25);


% Set up reaction parameters
obj.rxn_stoich(1,1,1) = -1; %species id, rxn id, subdomain id
obj.rxn_stoich(1,2,1) = 1; %species id, rxn id, subdomain id
obj.rxn_act_energy(1,1) = 25000; % rxn id, subdomain id
obj.rxn_act_energy(2,1) = 55000; % rxn id, subdomain id
obj.rxn_rate_const(1,1) = 1e2;
obj.rxn_rate_const(2,1) = 1e-2;
obj.rxn_powers(1,1,1) = 1;
obj.rxn_powers(2,2,1) = 1;
obj.rxn_enthalpy(1,1) = -0.5e6;
obj.rxn_enthalpy(2,1) = 0.5e7;

obj.rxn_stoich(2,1,1) = 1; %species id, rxn id, subdomain id
obj.rxn_stoich(2,2,1) = -1; %species id, rxn id, subdomain id
obj.mobile_spec_idx(2,1,1) = 0; % 2nd species is immobile 

% Call this before setting BCs, but after setting parameters
obj.set_coefficients();

% Setup input BCs
%       5,6,7,8 - Injection of stuff
%       9,10,11,12 - Pulling water out
inbound_set = [5,6,7,8, 9,10,11,12];
velocity_set = [0.005,0.005,0.005,0.0005, -0.005,-0.005,-0.005,-0.005];
temperature_set = [315,315,315,315, 298,298,298,298];
concentration_matrix = [0.6,0.6,0.6,0.6, 0.0,0.0,0.0,0.0;
                        0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0]; %NxID matrix of concentrations (2x8)
obj.set_input_boundaries(inbound_set, velocity_set, ...
                    temperature_set, concentration_matrix);


% Setup output BCs
outbound_set = [1,2,3,4];
pressure_set = [101350,101350,101350,101350];
obj.set_output_boundaries(outbound_set, pressure_set);



% Set initial conditions 
subdomain_set = [1];
pressure_set = [101350];
temperature_set = [298];
concentration_matrix = [0.01;
                        0.02];
obj.set_initial_conditions(subdomain_set,pressure_set,temperature_set, concentration_matrix);

% Setup simulation and run
nsteps=20;
t_span = linspace(0,200,nsteps);
results = obj.solve_system(t_span);
u = results.NodalSolution;


% Plot results
f1 = figure;
pdeplot(obj.model,"XYData", u(:,1,end),"ZData",u(:,1,end) ,Mesh="on", ColorMap="jet")
saveas(f1,'output/Gifs/pressure_test02_injection_well.png');
close(f1);

% Create gif for temperature
first_command = 'pdeplot(model,"XYData",u(:,2,INDEX),"ZData",u(:,2,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,2,:)));','umin = min(min(u(:,2,:)));',...
                'axis([-10 15 -15 10 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel T;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output2 = create_gif(command_set,variable_set,index_limit,1,'temp_solution_test02_injection_well');



% Create gif for concentration A
first_command = 'pdeplot(model,"XYData",u(:,3,INDEX),"ZData",u(:,3,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,3,:)));','umin = min(min(u(:,3,:)));',...
                'axis([-10 15 -15 10 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel C_A;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output3 = create_gif(command_set,variable_set,index_limit,1,'conc_solution_test02_CA_injection_well');


% Create gif for concentration B
first_command = 'pdeplot(model,"XYData",u(:,4,INDEX),"ZData",u(:,4,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,4,:)));','umin = min(min(u(:,4,:)));',...
                'axis([-10 15 -15 10 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel C_B;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output4 = create_gif(command_set,variable_set,index_limit,1,'conc_solution_test02_CB_injection_well');
