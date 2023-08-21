%% Test 4 - Multi-domain properties and reactions

% Create geometry for the pde
%   2D geometry 

% Rectangle is code 3, 4 sides (in m)
%          x1,x2,x3,x4       y1,y2,y3,y4
R1 = [3,4, 0, 0.05, 0.05, 0,        0, 0, 0.025, 0.025]';
R2 = [3,4, 0.05, 0.1, 0.1, 0.05,        0, 0, 0.025, 0.025]';
geom = [R1,R2];

% Names for the two geometric objects
ns = (char('R1','R2'))';

% Set formula
sf = 'R1 + R2';

% Create geometry
% E1 - bottom, E2 - right, E3 - top, E4 - left
g = decsg(geom,sf,ns);

% Plot results
f1 = figure;
pdegplot(g,"EdgeLabels","on","FaceLabels","on")
saveas(f1,'output/Gifs/pressure_test03_subdomains.png');
close(f1);

% BCs - boundary id 1 = input, positive velocity
%       boundary id 2 = output,
%       face id 1 = subdomain 1
%       face id 2 = subdomain 2

% Rxns
%       Subdomain 1
%           A --> B
%           B --> C
%       Subdomain 2
%           C --> D
%           D --> C
%
% Mobility: B and D are immobile, A and C are mobile 

% Create instance of PDE object 
%                   4 chemical species, 4 reactions, 2 subdomains
obj = porous_flow_2D(4, 4, 2);

% Set geometry based on edges and setup functions
obj.set_geometry_from_edges(g,"linear",0.005);

% Set up reaction parameters: r1 --> A & B in Sub01
obj.rxn_stoich(1,1,1) = -1; %species id, rxn id, subdomain id
obj.rxn_act_energy(1,1) = 25000; % rxn id, subdomain id
obj.rxn_rate_const(1,1) = 1e3;
obj.rxn_powers(1,1,1) = 1;
obj.rxn_enthalpy(1,1) = -5e7;
obj.rxn_stoich(2,1,1) = 1; %species id, rxn id, subdomain id

% Set up reaction parameters: r2 --> B & C in Sub01
obj.rxn_stoich(2,2,1) = -1; %species id, rxn id, subdomain id
obj.rxn_act_energy(2,1) = 55000; % rxn id, subdomain id
obj.rxn_rate_const(2,1) = 1e8;
obj.rxn_powers(2,2,1) = 1;
obj.rxn_enthalpy(2,1) = 3e3;
obj.rxn_stoich(3,2,1) = 1; %species id, rxn id, subdomain id

% Set up reaction parameters: r3 --> C & D in Sub02
obj.rxn_stoich(3,3,2) = -1; %species id, rxn id, subdomain id
obj.rxn_act_energy(3,2) = 15000; % rxn id, subdomain id
obj.rxn_rate_const(3,2) = 5e1;
obj.rxn_powers(3,3,2) = 1;
obj.rxn_enthalpy(3,2) = 2e7;
obj.rxn_stoich(4,3,2) = 1; %species id, rxn id, subdomain id

% Set up reaction parameters: r4 --> D & C in Sub02
obj.rxn_stoich(4,4,2) = -1; %species id, rxn id, subdomain id
obj.rxn_act_energy(4,2) = 75000; % rxn id, subdomain id
obj.rxn_rate_const(4,2) = 1e3;
obj.rxn_powers(4,4,2) = 1;
obj.rxn_enthalpy(4,2) = -5e3;
obj.rxn_stoich(3,4,2) = 1; %species id, rxn id, subdomain id

obj.mobile_spec_idx(2,1,1) = 0; % B species is immobile 
obj.mobile_spec_idx(2,1,2) = 0; % B species is immobile 
obj.mobile_spec_idx(4,1,1) = 0; % D species is immobile 
obj.mobile_spec_idx(4,1,2) = 0; % D species is immobile 


% Setup subdomain physical constants (Much more can be altered than shown
% here, but this is just for an example)
obj.bulk_porosity(1,1,1) = 0.33;
obj.bulk_porosity(1,1,2) = 0.45;
obj.particle_diameter(1,1,1) = 4e-3;
obj.particle_diameter(1,1,2) = 6e-3;
obj.char_len = obj.particle_diameter;
obj.bulk_solids_dens(1,1,1) = 2100;
obj.bulk_solids_dens(1,1,2) = 2500;
obj.bulk_solids_Cp(1,1,1) = 1.715e3;
obj.bulk_solids_Cp(1,1,2) = 1.115e3;


% Call this before setting BCs, but after setting parameters
obj.set_coefficients();

% Create BCs
% Setup input BCs
%       1
inbound_set = [1];
velocity_set = [0.0005];
temperature_set = [298];
concentration_matrix = [0.6;
                        0.0;
                        0.0;
                        0.0]; %NxID matrix of concentrations (4x1)
obj.set_input_boundaries(inbound_set, velocity_set, ...
                    temperature_set, concentration_matrix);

% Setup output BCs
outbound_set = [2];
pressure_set = [101350];
obj.set_output_boundaries(outbound_set, pressure_set);

% Set initial conditions (for each subdomain)
subdomain_set = [1,2];
pressure_set = [101350,101350];
temperature_set = [298,298];
concentration_matrix = [0.0, 0.0;
                        0.0, 0.0;
                        0.0, 0.0;
                        0.0, 0.0];  % NxID (4x2)
obj.set_initial_conditions(subdomain_set,pressure_set,temperature_set, concentration_matrix);

% Setup simulation and run
nsteps=80;
t_span = linspace(0,250,nsteps);
results = obj.solve_system(t_span);
u = results.NodalSolution;


% Plot results
f1 = figure;
pdeplot(obj.model,"XYData", u(:,1,end),"ZData",u(:,1,end) ,Mesh="on", ColorMap="jet")
saveas(f1,'output/Gifs/pressure_test03_multidomain.png');
close(f1);


% Create gif for temperature
first_command = 'pdeplot(model,"XYData",u(:,2,INDEX),"ZData",u(:,2,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,2,:)));','umin = min(min(u(:,2,:)));',...
                'axis([0 0.1 0 0.025 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel T;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output2 = create_gif(command_set,variable_set,index_limit,1,'temp_solution_test03_multidomain');



% Create gif for concentration A
first_command = 'pdeplot(model,"XYData",u(:,3,INDEX),"ZData",u(:,3,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,3,:)));','umin = min(min(u(:,3,:)));',...
                'axis([0 0.1 0 0.025 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel C_A;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output3 = create_gif(command_set,variable_set,index_limit,1,'conc_solution_test03_CA_multidomain');


% Create gif for concentration B
first_command = 'pdeplot(model,"XYData",u(:,4,INDEX),"ZData",u(:,4,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,4,:)));','umin = min(min(u(:,4,:)));',...
                'axis([0 0.1 0 0.025 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel C_B;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output4 = create_gif(command_set,variable_set,index_limit,1,'conc_solution_test03_CB_multidomain');



% Create gif for concentration C
first_command = 'pdeplot(model,"XYData",u(:,5,INDEX),"ZData",u(:,5,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,5,:)));','umin = min(min(u(:,5,:)));',...
                'axis([0 0.1 0 0.025 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel C_C;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output5 = create_gif(command_set,variable_set,index_limit,1,'conc_solution_test03_CC_multidomain');



% Create gif for concentration D
first_command = 'pdeplot(model,"XYData",u(:,6,INDEX),"ZData",u(:,6,INDEX),"ZStyle","continuous",ColorMap="jet");';
command_set = [first_command,'umax = max(max(u(:,6,:)));','umin = min(min(u(:,6,:)));',...
                'axis([0 0.1 0 0.025 umin umax]);','clim([umin umax]);',...
                'xlabel x;','ylabel y;','zlabel C_D;'];

variable_set = cell(2,2);
variable_set{1,1} = 'u';
variable_set{1,2} = u;
variable_set{2,1} = 'model';
variable_set{2,2} = obj.model;

index_limit = nsteps;
output6 = create_gif(command_set,variable_set,index_limit,1,'conc_solution_test03_CD_multidomain');


% Plot "breakthrough" concentration of A and C

% Interpolate solution on a line
%
%   line from point (0.1,0) to (0.1, 0.025)  [exit]
%                   (x1, y1) to (x2, y2)
xq = 0.1*ones(1,10);
yq = linspace(0,0.025,10);

CA_idx = 3;
CC_idx = 5;
CA_intrp = interpolateSolution(results,xq,yq, CA_idx, 1:length(t_span));
CC_intrp = interpolateSolution(results,xq,yq, CC_idx, 1:length(t_span));

% Interpolate all concentrations along x-axis
xq2 = linspace(0,0.1,50);
yq2 = 0.01*ones(1,50);
Cxset_intrp = interpolateSolution(results,xq2,yq2, 3:6, 1:length(t_span));

% CA_intrp will be a 10x1x80 (SpacePoints x Vars x TimePoints) structure 
breakthroughA = zeros(1,1,length(t_span));
breakthroughC = zeros(1,1,length(t_span));
for i=1:length(t_span)
    breakthroughA(1,1,i) = mean(CA_intrp(:,1,i));
    breakthroughC(1,1,i) = mean(CC_intrp(:,1,i));
end
times = results.SolutionTimes;
CA_out = breakthroughA(1,1,:);
CA_out = reshape(CA_out,[1,length(t_span)]);
CC_out = breakthroughC(1,1,:);
CC_out = reshape(CC_out,[1,length(t_span)]);

f2 = figure;
plot(times,CA_out, times,CC_out)
xlabel('Time (s)')
ylabel('Exit Concentration (mol/m^3)')
legend('C_A','C_C')
saveas(f2,'output/Gifs/breakthrough_conc_test03_multidomain.png');
close(f2);

f3 = figure;
plot(xq2,Cxset_intrp(:,:,end))
xlabel('Distance in x (m)')
ylabel('Concentration (mol/m^3)')
legend('C_A','C_B','C_C','C_D')
saveas(f3,'output/Gifs/steadystate_profile_conc_test03_multidomain.png');
close(f3);
